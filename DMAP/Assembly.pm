package DMAP::Assembly;
use PDL;
#use PDL::Fit::LM;
use PDL::LiteF;
use PDL::Fit::Levmar;
use PDL::NiceSlice;
use DMAP::Molecule;
use strict;
no strict 'refs';
=head DMAP::Assembly.pm

Wrapper object for a pseudomolecule analysis. This will match the model to a polynomial and calculate individual variations for the various component molecules.


=head2 new()

Create a new instance of a physical assembly. 

 my $ass=DMAP::Assembly->new(\%options);

options are:

=over

=item name (optional) name to use for the molecule.

=back

=cut


sub new {
    my ($class, $opts)=@_;
    my $self={
	name=>"",
	molecules=>{},
	ordered_molecules=>[],
	fitparams=>{},
	sorted=>0,
	markers=>{}, # hash to hold reference of molecule containing each marker for map updates..
	length=>0,
	minmappos=>1000,
	maxmappos=>0
    };
    if ($opts) {
	foreach my $k (keys %$opts){
	    if ($opts->{$k}){
		$self->{$k}=$opts->{$k};
	    }
	}
    }
    bless $self, $class;
    return $self;
}

=head2 addMolecule()

Add a molecule entry to the assembly. 

C<< $ass->addMolecule($start, $length, $name [,$orient, [$offset]] 
$start is the start position for the molecule, $length is the molecule length and $name is the molecule identifier. Optionally an orientation can be given. The molecule will be recorded as flipped if $orient is true. $offset is the length of the molecule before the start base in the molecule fragment (default 0).

=cut

sub addMolecule {
    my ($self,$start, $length, $name, $orient,$offset )=@_;
    if ($start,$length, $name) {
	my $o=$orient?1:0;
	my $offset=$offset?$offset:0;
	my $m=DMAP::Molecule->new({name=>$name, length=>$length, orient=>$o, offset=>$offset});
	if (ref($m) eq 'DMAP::Molecule'){
	    if (exists($self->{molecules}{$name})){
	    	my $frags=scalar @{$self->{molecules}{$name}};
	    	$m->fragment($frags);
		    push @{$self->{molecules}{$name}},{name=>$name, start=>$start, molecule=>$m, origstart=>$start, fragment=>$frags, origorient=>$o, offset=>$offset, length=>$length};
	    } else {
	    	@{$self->{molecules}{$name}}=({name=>$name, start=>$start, molecule=>$m, origstart=>$start, fragment=>0, origorient=>$o, offset=>$offset, length=>$length});
	    }
	    #print STDERR "adding $name at $start\n";
	    $self->{sorted}=0;
	    if (($start+$length -1) > $self->{length}){
		$self->{length}=$start+$length -1;
	    }
	}
    }    
}




=head2 _sortmol

provides a sorted list of molecules.

=cut

sub _sortmol {
    my ($self)=@_;
    #TODO change to use getAllMolecules instead of direct access. #DONE
    my @mol=();
    foreach my $k (values %{$self->{molecules}}){
    	push @mol, @{$k};
    }
    @{$self->{ordered_molecules}}=sort {$a->{start} <=>$b->{start}} @mol;
    $self->{sorted}=1;
}

=head2 addMarker()

Adds a link between the genetic and physical map. The position should be that calculated with reference to the physical map.

$ass->addMarker($name, $geneticmapposition, $physicalmapposition, $type);

=cut

sub addMarker {

    my ($self, $name, $mappos, $molpos, $type)=@_;

#binary search molecules to find one to put the marker in.
#then add marker, having tested the molecute to see if it is flipped.
    unless ($name && $mappos && $molpos && $type) {
	return;
    }
    unless (scalar @{$self->{ordered_molecules}}==scalar keys %{$self->{molecules}}){
	$self->_sortmol();
    }
    if ($molpos <0 || $molpos > $self->{length}) {
	warn "Marker position not on pseudomolecule\n";
	return;
    }
    if ($mappos > $self->{maxmappos}){
	$self->{maxmappos}=$mappos;
    }
    if ($mappos < $self->{minmappos}){
	$self->{minmappos}=$mappos;
    }
    my $molc=scalar keys %{$self->{molecules}};
    my $startpos=int($molc*$molpos/$self->{length});
    
    my $mn=$self->{ordered_molecules}[$startpos];
    #print STDERR "$mn, $startpos\n";
    my $lastmin=0; 
    my $lastmax=$molc-1;
    while ($molpos < $self->{molecules}{$mn}{start} || 
	   $molpos > $self->{molecules}{$mn}{start} +$self->{molecules}{$mn}{molecule}->length() 
	   && ! ($startpos==$lastmin || $startpos==$lastmax )) {
#	    print STDERR "checking start position\n";
#	    print STDERR "$name, $mappos, $molpos, ",$self->{molecules}{$mn}{start},",",$self->{molecules}{$mn}{molecule}->length(),",$mn\n";
	if ($molpos < $self->{molecules}{$mn}{start}){
	    $lastmax=$startpos;
	    $startpos=int($lastmin-($lastmin+$startpos)/2);
	}else{ 
	    $lastmin=$startpos;
	    $startpos=int(($startpos+$lastmax)/2);
	}
	$mn =$self->{ordered_molecules}[$startpos];
    }
    while ($molpos < $self->{molecules}{$mn}{start} ){
	$startpos++;
	$mn =$self->{ordered_molecules}[$startpos];
    }
    while ( $molpos > $self->{molecules}{$mn}{start} +$self->{molecules}{$mn}{molecule}->length()) {
 	$startpos--;
#	print STDERR "$mn, $startpos,$molpos, $self->{molecules}{$mn}{start},",$self->{molecules}{$mn}{molecule}->length(),",",$self->{ordered_molecules}[$startpos],"\n";
	$mn =$self->{ordered_molecules}[$startpos];
	exit unless ($self->{ordered_molecules}[$startpos]);
    }	
    my $mp=1+$molpos-$self->{molecules}{$mn}{start};
    if ($self->{molecules}{$mn}{molecule}->isFlipped()){
	$mp=1+$self->{molecules}{$mn}{length}-$mp;
    }
    my $mapname="default";
    $self->{markers}{$name}=$mn;
    $self->{molecules}{$mn}{molecule}->addMarker($name, $type,$mp, $mappos,$mapname ); # this adds to the default map
}

=head2 addGFFMarker()

Adds a link between the genetic and physical map. The position should be that calculated with reference to a molecule in the physical map.

$ass->addMarker($molname, $markername, $geneticmapposition, $physicalpositioninmolecule, $type);

=cut

sub addGFFMarker {

#    my ($self, $molname, $name, $mappos, $molpos, $type, $mapname)=@_;
    my ($self, $molname, $name, $molpos, $type)=@_;

#binary search molecules to find one to put the marker in.
#then add marker, having tested the molecute to see if it is flipped.
#    unless ($molname && $name && $mappos && $molpos&& $type) {
    unless ($molname && $name && $molpos && $type) {
		return;
    }
#    unless($mapname){
#    	$mapname="default";
#    }
#TODO change to refer to $self->{sorted} #DONE
    unless ($self->{sorted}){
	$self->_sortmol();
    }
    unless (exists($self->{molecules}{$molname})) {
		warn "Marker position not found in pseudomolecule\n";
		return;
    }

#    if ($mappos > $self->{maxmappos}){#
#		$self->{maxmappos}=$mappos;
#    }
#    if ($mappos < $self->{minmappos}){#
#		$self->{minmappos}=$mappos;
#    }
    $self->{markers}{$name}=$molname;
#    $self->{molecules}{$molname}{molecule}->addMarker($name, $type,$molpos, $mappos, $mapname);
    #TODO need to add to the appropriate molecule #DONE
    my $f=0;
    my @frags=sort {$a->{start}<=> $b->{start}} @{$self->{molecules}{$molname}};
    while ($f<scalar @frags && $molpos > $frags[$f]->{offset}+$frags[$f]->{length}){
    	print STDERR "FRAGS: ".join("::", $f, scalar @frags , $molpos , $frags[$f], $frags[$f]->{offset}, $frags[$f]->{length})."\n";
    	$f++;
    }
    if ($f < scalar @frags && $molpos >$frags[$f]->{offset}){
	    $frags[$f]->{molecule}->addMarker($name, $type,$molpos);
    }	else {
    	warn "Marker position not found in fragment\n";
    }
}

sub addMapPos {
	#TODO need to check references here #DONE
	my ($self, $map, $marker,$pos)=@_;
	my $success=0;
	if (exists($self->{markers}{$marker}) && 
		exists($self->{molecules}{$self->{markers}{$marker}})){
			foreach my $m (@{$self->{molecules}{$self->{markers}{$marker}}}){
				if ($m->{molecule}->hasMarker($marker)){
					$m->{molecule}->addMapPos($map, $marker, $pos);
					$success++;
					
				}
			}
	}
	return $success;
}

=head2 getMarkers

returns a list of all markers in the pseudomolecule.

my @marks=$ass->getMarkers();

C<<<foreach my $m (@marks){
    my $mappos=$m->{mappos};
    my $molpos=$m->{molpos};
    my $name=$m->{name};
    my $type=$m->{type};
}>>>

=cut

sub getMarkers {
    my ($self, $mapname)=@_;
    unless($mapname){
    	$mapname="default";
    }
    print STDERR "retreiving marker positions for map $mapname\n";
    my @markers=();
    #TODO use sorted molecule list #DONE
    unless ($self->{sorted}){
    	$self->_sortmol();
    }
    foreach my $mol (@{$self->{ordered_molecules}}){
    	print STDERR "Getting markers from $mol ".$mol->{name}."\n";
		foreach my $m ($mol->{molecule}->getMarkers($mapname)){
	    	my $molpos=$m->{molpos}+$mol->{start};
	    	push @markers, {molpos=>$molpos, avemappos=>$m->{mappos}, mappos=>$m->{mappos}, type=>$m->{type}, name=>$m->{name}};
		}
    }
    return @markers;
}

=head2 fit()

fit the data to the original model, returning that fit.

my %fit=$ass->fit($mapname);

=cut

sub fit {
    my ($self,$mapname)=@_;
    unless($mapname){
    	$mapname="default";
    }
    unless (exists($self->{originalfit})){
	%{$self->{originalfit}}=$self->_getFit($mapname);
    }
    my %fit=();
    foreach my $k (keys %{$self->{originalfit}}){
	$fit{$k}=$self->{originalfit}{$k};
    }
    return %fit;
}

=head2 bestfit

returns the optimised fit, calculating if necessary.

my %bestfit=$ass->bestfit(); 

=cut

sub bestfit {
    my ($self, $mapname)=@_;
    unless ($mapname) {$mapname="default";}
    unless (exists($self->{currentfit}{chi2})){
	$self->optimiseFlip($mapname);
    }
    my %fit=();
    foreach my $k (keys %{$self->{currentfit}}){
	$fit{$k}=$self->{currentfit}{$k};
	#print STDERR "retrieving $k $fit{$k} for Bestfit\n";
    }
    return %fit;
}

=head2 report

creates a hash containing a detaield report on fitting for all the molecules and the overall assembly.


=cut

sub report {
    my ($self, $mapname)=@_;
#Build a report for all the changes.
    unless ($mapname) { $mapname="default";}
    my %report=();
    %{$report{assembly}{original}}=$self->fit($mapname);
    %{$report{assembly}{final}}=$self->bestfit($mapname);
    @{$report{molecules}}=();
    #TODO use ordered_molecules DONE
    unless ($self->{sorted}){
    	$self->_sortmol();
    }
    foreach my $mr (@{$self->{ordered_molecules}}){
    	my $mol=$mr->{name};
		my %m=%{$mr};
		my $mstats={
	    name=>join("_",$mol, $m{fragment}),
	    start=>$m{start},
	    length=>$m{molecule}->length(),
	    offset=>$m{offset},
	    originalorient=>$m{originalorient}?"reverse":"forward",
	};
	if (exists($m{chi2})){
	    $mstats->{originalchi2}=$m{chi2};
	    $mstats->{originalvar}=$m{var};
	    $mstats->{neworient}=$m{molecule}->isFlipped()?"reverse":"forward";
	    my %newfit=$self->_fitMolecule($mr, $mapname);
	    $mstats->{newchi2}=$newfit{chi2};
	    $mstats->{newvar}=$newfit{var};
	    $mstats->{markers}=$newfit{markers};
	} else {
	    $mstats->{originalchi2}="";
	    $mstats->{originalvar}="";
	    $mstats->{neworient}=$mstats->{originalorient};
	    $mstats->{newchi2}="";
	    $mstats->{newvar}="";
	    $mstats->{markers}=0;
	}
	push @{$report{molecules}},$mstats;
    }
    return \%report;

}

=head2 optimiseFlip()

Fits the data, if not already done, then tries to improve the fit by sequentially flipping the molecules starting with the highest chi2

=cut

sub optimiseFlip {
    my ($self, $mapname)=@_;
    unless ($mapname){ $mapname="default";}
    my %fit=$self->fit($mapname);
    foreach my $k (keys %fit){
	$self->{currentfit}{$k}=$fit{$k};
    }
    my %molfits=();
    #TODO tricky as passes a molecule reference. should a hash reference be passed instead?
    foreach my $molref ( @{$self->{ordered_molecules}}){
		%fit=$self->_fitMolecule($molref, $mapname);
		if (exists($fit{markers})){
	    	foreach my $p (keys %fit){
				$molfits{join("_",$molref->{name}, $molref->{fragment})}{$p}=$fit{$p};
	    	}
	    #TODO need to find way for molecule to reference itself properly.
	    	$molref->{chi2}=$fit{chi2};
	    	$molref->{var}=$fit{var};
	    #print STDERR "FIT: $mol $fit{chi2} $fit{var}\n";
		}
    }
    my $changes=1;
    my $iter=1;
    while ($changes && $iter<10){ #iterate up to ten times to convergence
		$changes=0;
    # now have an initial fit and a list of molecules with their fits to the global fit.
		foreach my $k (keys %molfits){
	    #print STDERR "MOL $k, $molfits{$k}{chi2}, $molfits{$k}{var} $molfits{$k}{markers}\n";
		}
		sub _chisort {  my $x=($molfits{$b}{chi2} <=> $molfits{$a}{chi2});  return int($x) }

#print STDERR "sort $a($molfits{$a}{chi2}), $b($molfits{$b}{chi2}) returns ", ($molfits{$b}{chi2} <=> $molfits{$a}{chi2}), "\n";
 		foreach my $fmn (sort _chisort keys %molfits ){ #take largest chi2 first
#	foreach my $fm ( keys %molfits ){ #take largest chi2 first
			unless($fmn=~/^(.+)_(\d+)$/){
				warn "Bad fragement identifier\n";
			}
			my $fm=$1;
			my $fn=$2;
		
		    if ($molfits{$fmn}{markers}>1){ # no point in flipping single marker molecules.
	    	#TODO split molecule to get fragment name. #DONE
		    	my $fmol=$self->{molecules}{$fm}[$fn];
				$fmol->{molecule}->flip();
				my %newfit=$self->_fitMolecule($fmol, $mapname);
				if ($newfit{chi2}< $molfits{$fmn}{chi2}){
		   # print STDERR "Flipped $fm improves chi2 from $molfits{$fm}{chi2} to $newfit{chi2}\n";
		    		$changes++;
				} else {
		   # print STDERR "Flipped $fm no improvement in chi2 from $molfits{$fm}{chi2} to $newfit{chi2}\n";
		    		$fmol->{molecule}->flip();
				}
	    	}
		}
		my %newfit=$self->_getFit($mapname);
		foreach my $k (keys %newfit){
	    	$self->{currentfit}{$k}=$newfit{$k};
		}
#	$changes=0;
		$iter++;
    }
    
}

sub getMolecule {
    my ($self, $name, $fragment)=@_;
    unless ($fragment){
    	$fragment=0;
    }
    if (exists($self->{molecules}{$name})){
    	#TODO deal with fragments. #DONE
		my %mol=(name=>$self->{molecules}{$name}[$fragment]{name},start=>$self->{molecules}{$name}[$fragment]{start},markers=>$self->{molecules}{$name}[$fragment]{markers},molecule=>$self->{molecules}{$name}[$fragment]{molecule});
		return %mol;
    } else {
    	warn "No such molecule $name\n";
    }
}

=head2 getAllMolecules

returns an array of hash of all molecules.

=cut

sub getAllMolecules {
    my ($self, $mapname)=@_;
    my @mols=();
    unless ($mapname) {$mapname="default";}
    #TODO deal with fragments. use ordered_molecules #DONE
    unless ($self->{sorted}) {
    	$self->_sortmol();
    }
    foreach my $m (@{$self->{ordered_molecules}}){
		if (0){
    		print STDERR  "checking $m\n";
			warn "ordered Molecule:\n";
    		foreach  my $k (keys %$m){
    			warn "$k: $m->{$k}\n";
    		}		
    	}
		push @mols,{name=>$m->{molecule}->{name},start=>$m->{origstart},markers=>$m->{markers},molecule=>$m->{molecule}};	
    }
    return @mols;
}


=head2 _fitMolecule

estimate mean error and mean square of error for one molecule. 

%fit=$self->_fitMolecule($molname,$mapname);

If there are markers it returns a hash  {chi2, var, markers} where chi2 is the mean square error, var is the mean error, and markers is the marker count.

=cut
#need to deal with fragment references in this (molname, fragnumber)
#TODO take molecule reference rather than name #DONE
sub _fitMolecule {
    my ($self, $molref, $mapname)=@_;
    unless($mapname){
    	$mapname="default";
    }
    unless (exists($self->{molecules}{$molref->{name}})){
		warn "trying to fit non-existent molecule ".$molref->{name}."\n";
		return;
    }
    my %molfits=();
    my %mol=%{$molref};
    my @mark=$mol{molecule}->getMarkers($mapname);
    if (@mark){
	my @x=();
	my @y=();
	foreach my $m (@mark){
	    push @x, $m->{molpos}+$mol{start};
	    push @y, $m->{mappos};
	}
	#print STDERR "fitting molecule with ",join(",",@x),":",join(",", @y),":",$self->{currentfit}, "\n";
	my ($chi2, $var)=_chi2(\@x, \@y, $self->{currentfit});
	%molfits=(chi2=>$chi2, var=>$var, markers=>scalar @mark);
    }
    return %molfits;
}
=head2 _getFit()

$self->_getFit($mapname)

internal method for calculating the current fit returns a hash of fit parameters.

=cut



#TODO need to deal with fragments. #DONE
sub _getFit {

    my ($self, $mapname)=@_;
    unless($mapname){
    	$mapname="default";
    }
    my %fit=();
    my @x=();
    my @y=();
    foreach my $mol (@{$self->{ordered_molecules}}){
	#print STDERR "Molecule $mol\n"; 
	foreach my $m (	$mol->{molecule}->getMarkers($mapname)){
	    #print STDERR "x,y: ",$m->{mappos}, ",",$self->{molecules}{$mol}{start}+$m->{molpos},"\n";
	    push @y, $m->{mappos};
	    push @x, $mol->{start}+$m->{molpos};

	}
    }
    unless (@x && @y) { warn "Data error - no data for this map $mapname\n"; return;}
    my $pdlx= pdl @x;

    my $pdly=pdl @y;
    #set preliminary parameter estimates
    #my $pdlf=pdl $self->{minmappos},exp((log10($self->{length})-2)*log(10)),-exp(2*(log10($self->{length})-2)*log(10)),exp(3*(log10($self->{length})-2)*log(10));
    my $pdlf=pdl 50.0,1.0,0.5,9.0,float($self->{length}); 
    my $fix= pdl 0,0,0,0,1;
    print STDERR "Fitting $pdlx, $pdly\n";
    #perform linear model fit.
    
    my $lf=levmar(P=>$pdlf,  X=>$pdly, T=>$pdlx, FUNC => sub {
	my ($p, $x, $t)=@_;
	my ($p0, $p1, $p2, $p3, $p4)=list $p;
	$x.=$p0 +$p1*sinh(($t - $p2*$self->{length})*$p3/$self->{length});
		  } 
	);   

    my @fits=list $lf->{P};
#    print levmar_report($lf);
#    print "output $lf->{P} @fits :: ".scalar @fits."\n";
    $fit{a0}=$fits[0];
    $fit{a1}=$fits[1];
    $fit{a2}=$fits[2];
    $fit{a3}=$fits[3];
    $fit{size}=$self->{length};
    $fit{report}=$lf;
    $fit{iters}=$lf->{ITS};
    ($fit{chi2}, $fit{var})=_chi2(\@x, \@y, \%fit);
#calculate chi squared for all points then per molecule.
    return %fit;

}

=head2 _chi2

Calculates Chi2 for @x and @y arrays fitted to the polynomial described by %fit

my $chi2=chi2(\@x, \@y, \%fit);

=cut

sub _chi2 {

    my ($xref, $yref, $fitref)=@_;
    my @var=();
    my $total2=0;
    my $total=0;
    #print STDERR join ",",keys %$fitref,"\n";
    unless (scalar @$xref == scalar @$yref) {
	warn "unequal length arrays\n";
	return;
    }
    my $i=0;
    for ( $i=0; $i < scalar @$xref; $i++){
#	print STDERR join ":", $i, $xref->[$i], $yref->[$i],$fitref->{a0},$fitref->{a1},$fitref->{a2},$fitref->{a3},"\n";
	my $v=$fitref->{a0} +$fitref->{a1}*sinh( $fitref->{a3}*($xref->[$i] - $fitref->{a2}*$fitref->{size})/$fitref->{size});
#	my $v=$yref->[$i]-($fitref->{a0}+ $fitref->{a1}*$xref->[$i]+ $fitref->{a2}*($xref->[$i]**2)+ $fitref->{a3}*($xref->[$i]**3));
	$total += $v;
	$total2 += ($v**2);
	$var[$i]=$v;
    }
    
    return $total2/$i, $total/$i;

}



=head2 _chrfit()

Internal method for passing to PDL::LMfit

=cut
 
sub _chrfit {
    # copied directly from the LMfit perldoc with minor changes
	   # leave this line as is
	   my ($x,$par,$ym,$dyda) = @_;

	   # $m and $b are fit parameters, internal to this function
	   # call them whatever make sense to you, but replace (0..1)
	   # with (0..x) where x is equal to your number of fit parameters
	   # minus 1
	   my ($a3,$a2, $a1, $a0) = map { $par->slice("($_)") } (0..3);

	   # Write function with dependent variable $ym,
	   # independent variable $x, and fit parameters as specified above.
	   # Use the .= (dot equals) assignment operator to express the equality 
	   # (not just a plain equals)
	   $ym .= $a3 * $x**3 + $a2 * $x**2 +$a1 * $x + $a0;

	   # Edit only the (0..1) part to (0..x) as above
	   my (@dy) = map {$dyda -> slice(",($_)") } (0..3);

	   # Partial derivative of the function with respect to first 
	   # fit parameter ($m in this case). Again, note .= assignment 
	   # operator (not just "equals")
	   $dy[0] .= $x**3;

	   # Partial derivative of the function with respect to next 
	   # fit parameter ($b in this case)
	   $dy[1] .= $x**2;

	   # Add $dy[ ] .= () lines as necessary to supply 
	   # partial derivatives for all floating paramters.
	   $dy[2] .= $x;
	   $dy[3] .= 1;


   }


1;
