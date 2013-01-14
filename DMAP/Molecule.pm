package DMAP::Molecule;

use strict;
=pod

=head B<DMAP::Molecule>

This holds information relating to a physical sequence scaffold. 

Parameters:

=over 1

=item B<name> the name of the molecule

=item B<length> Molecule length

=item B<orient> the orientation of the molecule. 0 forward, 1 reversed. 

=back

Molecule holds an internal array of markers which contain name, type, molecule and map positions.

=head2 new

Create a new DMAP::Molecule instance

C<<< my $mol=DMAP::Molecule->new({length=>1234, 
    name=>"My Molecule", 
    orient=>0, 
    markers=>[]}); >>>

=cut


sub new {
    my ($class,$par)=@_;
    my $self={markers=>{},
	      orient=>0, # 0 for forwards 1 for reverse
	      name=>"",
	      offset=>0, #if this is a fragment of the molecule, the start coordinate of the fragment
	      fragment=>0,	
	      length=>0
    };

    if ($par->{name}){
		$self->{name}=$par->{name};
    }
    if ($par->{offset}){
    	$self->{offset}=$par->{offset};
    }
    if ($par->{fragment}){
    	$self->{fragment}=$par->{fragment};
    }
    if ($par->{length}){
		$self->{length}=$par->{length};
    } else {
		warn("trying to create a zero length molecule\n");
		return;
    }
    
    if ($par->{markers}){
	foreach my $m (@{$par->{markers}}){
	    $self->{markers}{$m->{name}}=$m;
	}
    }
    bless $self, $class;
    return $self;
}

=head2 addMarker()

Add a marker instance to the DMAP::Molecule.

C<< $mol->addMarker($name, $type, $molpos [, $mappos[ ,$mapname]] ); >>

=cut

sub addMarker {
    my ($self, $name, $type, $molpos, $mappos, $mapname)=@_;
    unless ($name && $type && $molpos ){
	warn("incomplete marker specification - not added\n");
    }
    if (! exists($self->{markers}{$name})){
		$self->{markers}{$name}= {name=>$name, type=>$type, molpos=>$molpos, mappos=>{}};
    }
    if ($mappos){
		unless($mapname){
	  	  $mapname="default";
		}
	$self->addMapPos($mapname, $name, $mappos);
    }
}

=head2 hasMarker()

boolean checks that a marker exists in the molecule

C<< $mol->hasMarker($markername); >>

=cut
sub hasMarker{
	my ($self, $name)=@_;
	if (exists($self->{markers}{$name})){
		return 1;
	}else{
		return 0;
	}
}


=head2 addMapPos()

Add another map position to the marker

C<< $mol->addMapPos($mapname, $markername, $mappos); >>

=cut

sub addMapPos {
    my ($self, $mapname, $markername, $mappos)=@_;
    if ($mapname && $markername && $mappos ) {
		if (exists($self->{markers}{$markername})){
		    $self->{markers}{$markername}{mappos}{$mapname}=$mappos;			
		} else {
		    warn "No such marker $markername in this molecule\n";
		}
    }
}



=head2 aveMapPos()

returns the average map position for the markers in this molecule.

my $pos=$mol->aveMapPos();

=cut

sub aveMapPos {
    my ($self, $mapname)=@_;
    unless ($mapname){
    	$mapname="default";
    }
    my $total=0;
    my $count=0;
    foreach my $m (values %{$self->{markers}}){
    	if (exists($m->{mappos}{$mapname})){
			$total+=$m->{mappos}{$mapname};	
			$count++;
    	}	
    }
    if ($count){
	return $total/$count;
    }
}

=head2 flip()

Invert the orientation of the molecule. When the molecule is inverted, the marker X positions (physical coordinates) are calculated on the reverse strand.

A subsequent call to C<flip()> will perform a subsequent inversion.

C<< $mol->flip(); >>


=cut

sub flip {
    my ($self)=@_;
    $self->{orient}=($self->{orient}+1)%2;
}

=head2 unflip()

reset the flip to forward strand.

C<< $mol->unflip(); >>

=cut

sub unflip {
    my ($self)=@_;
    $self->{orient}=0;
}

=head2 fragment()
 Can be used to set or get the fragment number
 
 <C $fragno=$mol->fragment(); >>
 
 <C $mol->fragment($fragno); >>
 
=cut

sub fragment {
	my ($self, $frag)=@_;
	if (defined($frag)){
		$self->{fragment}=$frag;
	}
	return $self->{fragment};
}

=head2 isflipped()

returns the flip status. 0 for forward, 1 for reverse.

C<< if ($mol->isflipped()) { ... } >>

=cut

sub isFlipped {
    my ($self)=@_;
    return $self->{orient};
}

=head2 getMarkers()

returns a list of markers as an array of anonymous hashes with molecule position calculated according to the orientation of the molecule. The hash fields are molpos, mappos, name and type.

 my @markers=$mol->getMarkers();
 foreach my $m (@markers){
    $x=$m->{molpos};
    $y=$m->{mappos};
    $name=$m->{name};
    $type=$m->{type};
 }

=cut

sub getMarkers {
    my ($self, $mapname)=@_;
    my @ma=();
    unless ($mapname){
    	$mapname="default";
    }
    foreach my $m (values %{$self->{markers}}){
    	if (exists($m->{mappos}{$mapname})) {	
			my $molpos=$self->{orient}?
				1+$self->{length}+$self->{offset}-$m->{molpos}:
				$m->{molpos}-$self->{offset};
			push @ma, {molpos=>$molpos, mappos=>$m->{mappos}{$mapname}, name=>$m->{name},type=>$m->{type}};
    	}
    }		
    return @ma;
}


sub length {
    my ($self)=@_;
    return $self->{'length'}-$self->{'offset'};

}



1;
