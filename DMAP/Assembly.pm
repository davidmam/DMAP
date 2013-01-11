package DMAP::Assembly;
use PDL;
use PDL::Fit::LM;
use DMAP::Molecule;
use strict;

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

C<< $ass->addMolecule($start, $length, $name [,$orient, [$offset]] );>>

$start is the start position for the molecule, $length is the molecule length and $name is the molecule identifier. Optionally an orientation can be given. The molecule will be recorded as flipped if $orient is true.

=cut

sub addMolecule {
    my ($self,$start, $length, $name, $orient,$offset )=@_;
    if ($start,$length, $name) {
	my $o=$orient?1:0;
	my $offset=$offset?$offset:0;
	my $m=DMAP::Molecule->new({name=>$name, length=>$length, orient=>$o, offset=>$offset});
	if (ref($m) eq 'DMAP::Molecule'){
	    $self->{molecules}{$name}={start=>$start, molecule=>$m, origstart=>$start, origorient=>$o};
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
    @{$self->{ordered_molecules}}=sort {$self->{molecules}{$a}{start} <=>$self->{molecules}{$b}{start}} keys %{$self->{molecules}};
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
    unless (scalar @{$self->{ordered_molecules}}==scalar keys %{$self->{molecules}}){
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
    $self->{molecules}{$molname}{molecule}->addMarker($name, $type,$molpos);
}

sub addMapPos {
	my ($self, $map, $marker,$pos)=@_;
	if (exists($self->{markers}{$marker}) && 
		exists($self->{molecules}{$self->{markers}{$marker}})){
			$self->{molecules}{$self->{markers}{$marker}}{molecule}->addMapPos($map, $marker, $pos);
	}
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
    foreach my $mol (sort {$a->{start}<=>$b->{start}} values %{$self->{molecules}}){
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
    foreach my $mol (sort {$self->{molecules}{$a}{start} <=> $self->{molecules}{$b}{start}} keys %{$self->{molecules}}){
	my %m=%{$self->{molecules}{$mol}};#
	my $mstats={
	    name=>$mol,
	    start=>$m{start},
	    length=>$m{molecule}->length(),
	    originalorient=>$m{originalorient}?"reverse":"forward",
	};
	if (exists($m{chi2})){
	    $mstats->{originalchi2}=$m{chi2};
	    $mstats->{originalvar}=$m{var};
	    $mstats->{neworient}=$m{molecule}->isFlipped()?"reverse":"forward";
	    my %newfit=$self->_fitMolecule($mol, $mapname);
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
    foreach my $mol (keys %{$self->{molecules}}){
	%fit=$self->_fitMolecule($mol, $mapname);
	if (exists($fit{markers})){
	    foreach my $p (keys %fit){
		$molfits{$mol}{$p}=$fit{$p};
	    }
	    $self->{molecules}{$mol}{chi2}=$fit{chi2};
	    $self->{molecules}{$mol}{var}=$fit{var};
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
 	foreach my $fm (sort _chisort keys %molfits ){ #take largest chi2 first
#	foreach my $fm ( keys %molfits ){ #take largest chi2 first
	    if ($molfits{$fm}{markers}>1){ # no point in flipping single marker molecules.
		$self->{molecules}{$fm}{molecule}->flip();
		my %newfit=$self->_fitMolecule($fm, $mapname);
		if ($newfit{chi2}< $molfits{$fm}{chi2}){
		   # print STDERR "Flipped $fm improves chi2 from $molfits{$fm}{chi2} to $newfit{chi2}\n";
		    $changes++;
		} else {
		   # print STDERR "Flipped $fm no improvement in chi2 from $molfits{$fm}{chi2} to $newfit{chi2}\n";
		    $self->{molecules}{$fm}{molecule}->flip();
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
    my ($self, $name)=@_;
    if (exists($self->{molecules}{$name})){
	my %mol=(name=>$self->{molecules}{$name}{name},start=>$self->{molecules}{$name}{start},markers=>$self->{molecules}{$name}{markers},molecule=>$self->{molecules}{$name}{molecule});
	return %mol;
    }
}

=head2 getAllMolecules

returns an array of hash of all molecules.

=cut

sub getAllMolecules {
    my ($self, $mapname)=@_;
    my @mols=();
    unless ($mapname) {$mapname="default";}
    foreach my $m (sort {$self->{molecules}{$a}{start} <=> $self->{molecules}{$b}{start}  } keys %{$self->{molecules}}) {
	print STDERR  "checking $m\n";
	
	push @mols,{name=>$m,start=>$self->{molecules}{$m}{start},markers=>$self->{molecules}{$m}{markers},molecule=>$self->{molecules}{$m}{molecule}};
	
    }
    return @mols;
}


=head2 _fitMolecule

estimate mean error and mean square of error for one molecule. 

%fit=$self->_fitMolecule($molname,$mapname);

If there are markers it returns a hash  {chi2, var, markers} where chi2 is the mean square error, var is the mean error, and markers is the marker count.

=cut

sub _fitMolecule {
    my ($self, $molname, $mapname)=@_;
    unless($mapname){
    	$mapname="default";
    }
    unless (exists($self->{molecules}{$molname})){
	warn "trying to fit non-existent molecule $molname\n";
	return;
    }
    my %molfits=();
    my %mol=%{$self->{molecules}{$molname}};
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




sub _getFit {

    my ($self, $mapname)=@_;
    unless($mapname){
    	$mapname="default";
    }
    my %fit=();
    my @x=();
    my @y=();
    foreach my $mol (keys %{$self->{molecules}}){
	#print STDERR "Molecule $mol\n"; 
	foreach my $m (	$self->{molecules}{$mol}{molecule}->getMarkers($mapname)){
	    #print STDERR "x,y: ",$m->{mappos}, ",",$self->{molecules}{$mol}{start}+$m->{molpos},"\n";
	    push @y, $m->{mappos};
	    push @x, $self->{molecules}{$mol}{start}+$m->{molpos};

	}
    }
    unless (@x && @y) { warn "Data error - no data for this map $mapname\n"; return;}
    my $pdlx= pdl @x;

    my $pdly=pdl @y;

    my $pdlf=pdl $self->{minmappos},exp((log10($self->{length})-2)*log(10)),-exp(2*(log10($self->{length})-2)*log(10)),exp(3*(log10($self->{length})-2)*log(10));
    print STDERR "Fitting $pdlx, $pdly\n";
    my ($yf, $pf, $cf, $if) =lmfit $pdlx, $pdly, 1, \&_chrfit, $pdlf;
    print STDERR "Fitted\n";
    $fit{yfit}=$yf;
    print STDERR $pf,"\n";
    $fit{a3}=$pf->index(0);
    $fit{a2}=$pf->index(1);
    $fit{a1}=$pf->index(2);
    $fit{a0}=$pf->index(3);
    $fit{cf}=$cf;
    $fit{iters}=$if;
   #$fit{mid}=-2*$fit{a2}/6*$fit{a3};
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
	my $v=$yref->[$i]-($fitref->{a0}+ $fitref->{a1}*$xref->[$i]+ $fitref->{a2}*($xref->[$i]**2)+ $fitref->{a3}*($xref->[$i]**3));
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
