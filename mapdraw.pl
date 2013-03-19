#!/usr/bin/perl -w

use PDF::API2;
use PDF::Table;
use Math::Trig;
use lib ".";
use DMAP::Assembly;
use strict;
=pod
=head mapdraw.pl

draws a genetic map aligned to a physical map.

=head2 file format

  MARKER,marker_name,marker_position,pseudomolecule_position,markertype
  MOLECULE,molecule_name,start,length
  LINK,type,colourname,linetype(dash|solid),font

options:
 -infile input file
 -outfile output file
 -height figure height in mm
 -width figure width in mm
-linewidth width of drawn lines in points
 -colour (multi) name:red,green,blue values 0-255
 -chrom chromosome colour
 -mol molecule colour
 -minfont minimum font size for markers (default 4)
 -crowd try to fit more marker labels in (default 2.5 - higher is more labels)
 -plotallmap plot every genetic map pos rather than binning.
 -nocount Don't include marker count.
 -genbin resolution of genetic map bin (default 0.1cM)
 -nosize don't add molecule size to the plot title.
 -seqbin resolution of sequence marker bin (only one marker of each type will be shown per bin) Default 100.
 -title plot title. 
 -unanchored Plot unanchored markers as grey bars on genetic map
 -trim Trim genetic map to maximum map position (otherwise 100 or max pos if greater). If -unanchored is not set then it trims to anchored markers only. 

Colours are X11/HTML/CSS colours.
Font is Helvetica or Times with options Roman Bold or Italic e.g. Helvetica Roman

=cut
#page definitions/resolution

# PDF::API2 native resolution is in points. This can be redefined into internal units.



use constant RES => 72/300;
use constant mm => 25.4 / 72;
use constant in => 1 / 72;
use constant pt => 1;



use Getopt::Long;
my $filename="drawing.pdf";
my @usercolours=();
my $chrcolour="default";
my $molcolour="default";
my $labelfont="Helvetica";
my $infile="";
my $height=297; #height in mm;
my $width=210; #width in mm
my $chrwidth=10; # chromosome width in mm - min 4
my $molwidth=6; # molecule width in mm - min 4
my $binsize=1.3;
my $chrheight=90; # percent of page height for a chromosome.
my $molfontsize=10;
my $markerfontsize=6;
my $minfont=4;
my $linewidth=0.5;
my $invert=0; # whether to have 0 at the bottom or top(default)
my $lineanglespace=4;
my $markerlabelwidth=32;
my $plotallmap=0;
my $crowd=2.5;
my $plotstart=0; # limits for physical molecule range to draw (0 no limit)
my $plotend=0;# limits for physical molecule range to draw (0 no limit)
my $maplabelwidth=15;
my $mollabelsize=55;
my $nosize=0;
my $molposbin=0;
my $plottitle="Chromosome Map";
my $genbin="0.1";
my $nocount=0;
my $debug=0;
my $correlate=""; # draw a correlation between the genetic and physical maps to the given file.
my $agpfile=""; # AGP file containing molecule definition.
my $gfffile=""; #GFF file containing marker positions on the molecules.
my $mapfile=""; #file containing marker positions on the molecules. Format "markername position .*"
my $mapname=""; #name of map (for versioning etc.)
my $chrtrim=0; # trim chromosome to map limits
my $plotunanchored=0; # plot unanchored linkage markers.

GetOptions(
    "infile=s"=>\$infile,
    "agpfile=s"=>\$agpfile,
    "gfffile=s"=>\$gfffile,
    "mapfile=s"=>\$mapfile,
    "mapname=s"=>\$mapname,
    "outfile=s"=>\$filename,
    "chrheight=i"=>\$chrheight,
    "height=i"=>\$height,
    "width=i"=>\$width,
    "invert"=>\$invert,
    "plotstart=i"=>\$plotstart,
    "plotend=i"=>\$plotend,
    "nosize"=>\$nosize,
    "plotallmap"=>\$plotallmap,
    "crowd=f"=>\$crowd,
    "seqbin=i"=>\$molposbin,
    "title=s"=>\$plottitle,
    "chrwidth=i"=>\$chrwidth,
    "colour=s"=>\@usercolours,
    "chrom=s"=>\$chrcolour,
    "correlation=s"=>\$correlate,
    "genbin=s"=>\$genbin,
    "minfont=i"=>\$minfont,
    "nocount"=>\$nocount,
    "linewidth=f"=>\$linewidth,
    "mol=s"=>\$molcolour,
    "unanchored"=>\$plotunanchored,
    "trim"=>\$chrtrim
    );

my $mapdp=0;
if (index($genbin,".")>=0){
    $mapdp=length($genbin)-index($genbin,".")-1;
}
#print STDERR join("\t",$mapdp, length($genbin), index($genbin,".")),"\n";

if ($molcolour eq 'default'){ $molcolour='black';}
if ($chrcolour eq 'default'){ $chrcolour='black';}
unless (-e $infile || ($agpfile && -e $agpfile && $mapfile && -e $mapfile && $gfffile && -e $gfffile)) {
    die "cannot find input file $infile\n";
}

my @markers=(); # list of markers to plot
my %colours=();
my %links=();
my @molecules=();
my $maxchrpos=0; # maximum genetic map distance
my $maxmolpos=0; # maximum physical map distance
my $minchrpos=1000;

my $assembly=DMAP::Assembly->new({name=>$plottitle});
my $mapfound=0;
my $agpfound=0;
my $gfffound=0;
my @unanchored=(); # holds unanchored genetic markers.

if ($agpfile && -e $agpfile) {
	open AGPFILE, "$agpfile" or die "Could not open AGP file: $!\n";
	my $mollen=0;
	while (my $agpline=<AGPFILE>){
		chomp $agpline;
		my @F=split /\t/,$agpline;
		next unless scalar @F >5;
		if ($F[4]=~m/[UN]/){
			$mollen+=$F[5];
		} else {
			my $wc=$F[8]eq'-'?1:0;
			my $fraglen=1+$F[7]-$F[6];

			$assembly->addMolecule($mollen+1,$fraglen,$F[5],$wc,$F[6]-1);
			print STDERR "adding molecule $F[7] $F[5]\n";
			$mollen+=$fraglen;


		}
	}
	close AGPFILE;
	$agpfound= $mollen;
	$maxmolpos=$mollen;
}

if ($gfffile && -e $gfffile) {
	open GFFFILE, "$gfffile" or die "Could not open GFF file: $!\n";
	my $gmap="default";
	if ($mapname && ! $mapfile) {
		$gmap=$mapname;
	}
	while (my $gff=<GFFFILE>){
		next if $gff=~m/^!/;
		chomp $gff;
		my @F=split /\t/, $gff;
		my %ext=();
		unless ($F[8]){
			$F[8]=$F[7];
		}
		if ($F[8]){
			foreach my $n (split (/; */ , $F[8])){
				my ($key, $value)= split( /=/,$n);
				$ext{$key}=$value;
			}
		}else{
			warn "Error with GFF line: $gff\n";
		}
		my %notes=();
		if (exists($ext{"Note"})){
		    #set default $notes{pos}
		    #$notes{pos}=-1;
		    foreach my $n (split(/, /, $ext{'Note'})){
			my ($key,$value)= split(/: /, $n);
			$notes{$key}=$value;
		    }
		}
		my $mmap=$gmap;
		#if ($notes{pos}==-1){ $mmap="fake";}
		next unless (exists($ext{ID}) #and exists($notes{pos}) 
			     and exists($notes{type}));
#		$assembly->addGFFMarker($F[0], $ext{ID},$notes{pos}, $F[3],$notes{type}, $mmap);
		$assembly->addGFFMarker($F[0], $ext{ID}, $F[3],$notes{type});

		if (exists($notes{pos}) && ! $mapfile ){
		    $assembly->addMapPos($mmap,$ext{ID}, $notes{pos});
		}
#		 $markername, $geneticmapposition, $physicalpositioninmolecule, $type);		
	}
	$gfffound=1;
	close GFFFILE;
}

if ($mapfile && -e $mapfile) {

    my $gmap=$mapname?$mapname:"default";
    open (MAPFILE, "$mapfile") or die "Cannto open map file $mapfile: $!\n";
    while (my $mapline=<MAPFILE>){
	next if $mapline=~m/^ *$/; #skip blank lines and comments.
	next if $mapline=~m/^;/; #skip blank lines and comments.
	next unless $mapline=~m/;/;
	#split lines and add to assembly
	my ($markername, $markerpos, $junk)= split /\s+/, $mapline, 3;
	if ($assembly->addMapPos($gmap, $markername, $markerpos)){
	    if ($markerpos > $maxchrpos){$maxchrpos=$markerpos;}
	    if ($markerpos < $minchrpos){$minchrpos=$markerpos;}
	    print STDERR "ADDING MAP POS $markername $markerpos $maxchrpos $minchrpos\n";
	} else {
	    if ($plotunanchored){
		push @unanchored, {"name"=>$markername, "pos"=>$markerpos};
		if ($markerpos > $maxchrpos){$maxchrpos=$markerpos;}
		if ($markerpos < $minchrpos){$minchrpos=$markerpos;}
	    }
	}
	
	
    }
    close MAPFILE;
    $mapfound=1;
}

if ($mapfound && $agpfound && $gfffound){  
    
    print STDERR "MAXCHR: $maxchrpos\n";

# need to build marker and molecule lists.
    #molecules
    foreach my $m ($assembly->getAllMolecules()) {
	#if (exists($m->{molecule}->{markers}{$m->{name}}{mappos}{$mapname})){
	    push @molecules, {name=>$m->{name}, start=>$m->{start}, length=>$m->{molecule}->{length} };
	    #print STDERR "name=>".$m->{name}.", start=>".$m->{start}.", length ".$m->{molecule}->{length}."\n";
	#} else {
	 #   warn "Marker ".$m->{name}." not in map $mapname\n";
	#}
    }
    #markers
    unless ($mapname){ $mapname="default";}
    @markers=$assembly->getMarkers($mapname);
    #push @markers, {name=>$f[0], avemappos=>avemappos($f[1]),mappos=>binmappos($f[1]), molpos=>$f[2], type=>$f[3]}; 

}
open (INFILE, $infile) or die "Error opening input file: $!\n";
while (my $line=<INFILE>){
    chomp $line;
    my @f=split /,/, $line;
    #print STDERR join(":",@f),"\n";
    my $type=shift @f;
	
    foreach ($type)  {
	/^MOLECULE$/ and do {  
	    push @molecules, {name=>$f[0], start=>$f[1], length=>$f[2] }; 
	    if ($f[1]+$f[2] > $maxmolpos) {$maxmolpos=$f[1]+$f[2];} 
	    $assembly->addMolecule($f[1], $f[2], $f[0]); 
	    last;};	
	/^LINK$/ and do {  $links{$f[0]}={ colour=>$f[1], line=>$f[2], font=>$f[3]}; last;};	
	/^MARKER$/ and do {  
	    push @markers, {name=>$f[0], avemappos=>avemappos($f[1]),mappos=>binmappos($f[1]), molpos=>$f[2], type=>$f[3]}; 
	    my $mmax=binmappos($f[1]); my $mmin=binmappos($f[1]);  
	    if ($mmax=~/(-?\d+)-(-?\d+)/){$mmax=$2; $mmin=$1} 
	    if ($mmax>$maxchrpos){$maxchrpos=$mmax;} 
		if ($mmin<$minchrpos){$minchrpos=$mmin;}  
	    $assembly->addMarker($f[0],avemappos($f[1]), $f[2],$f[3]);   
	    last ;};	      
    }
}

close INFILE;


my $report=$assembly->report($mapname);

foreach my $col (@usercolours){
    my ($n, $r, $g, $b)=split/[:,]/, $col;
    if ($n && $r && $b && $g) {
	$colours{$n}={red=>$r, green=>$g, blue=>$b};
    }
}


my $maxheight =$chrheight*$height/100;
my $basemargin=$height*(100-$chrheight)/200;

my $chrx=($width/5)-($chrwidth/2);
my $molx=(2*$width/5)-($molwidth/2);
my $mollabelwidth=$width-$mollabelsize;
my $pdf=PDF::API2->new( -file =>$filename);

die "could not create PDF\n" unless $pdf;

# set up font definitions.
my %font = (
    Helvetica => {
	Bold   => $pdf->corefont( 'Helvetica-Bold',    -encoding => 'latin1' ),
	Roman  => $pdf->corefont( 'Helvetica',         -encoding => 'latin1' ),
	Italic => $pdf->corefont( 'Helvetica-Oblique', -encoding => 'latin1' ),
    },
    Times => {
	Bold   => $pdf->corefont( 'Times-Bold',   -encoding => 'latin1' ),
	Roman  => $pdf->corefont( 'Times',        -encoding => 'latin1' ),
	Italic => $pdf->corefont( 'Times-Italic', -encoding => 'latin1' ),
    }
    );
my $moltextfont=$font{'Helvetica'}{'Bold'};
my $titletextfont=$font{'Helvetica'}{'Bold'};
my $titlefontsize=14;
my $titlecolour='black';
# define colours here.


#get page for main figure

my $page=$pdf->page;
$page->mediabox(int($width/mm), int($height/mm));

#get page for correlation plot.

my $cpage=$pdf->page;
$cpage->mediabox(int($width/mm), int($height/mm));

# get page for correlation plot with modelled curve.

my $dpage=$pdf->page;
$dpage->mediabox(int($width/mm), int($height/mm));

# get page for tabulated reports.

my $rpage=$pdf->page;
$rpage->mediabox(int($width/mm), int($height/mm));



my $caxiscolour='black';
my $cbasecolour='darkgray';
# draw a box and title. Plot will be square
#plotsize - box should be about 70%  of total available width


my $cplotsize=$width * 0.7;
my $cplotyaxis=$width * 0.2;
my $cplotxaxis=($height-$cplotsize)/2;
my $cplotystart=0<$minchrpos?0:$minchrpos;

my $cplotyscale=($maxchrpos>100?$maxchrpos:100)-$cplotystart;
if ($chrtrim){
    $cplotystart=$minchrpos;
    $cplotyscale=$maxchrpos-$cplotystart;
}
# plot the figure outline for the correlation and modelled plots.
plotfig($cpage);
plotfig($dpage);

sub plotfig { #plots the figure outline on the given page.
    my ($cpage)=@_;
    my $cplot=$cpage->gfx;
    if ($linewidth){
	my $lw=$cplot->linewidth($linewidth);
    
	#foreach my $k (keys %$lw){ 
	#    print STDERR "LINEWIDTH: $k $lw->{$k}\n";
	#}
    }
    my $ctext=$cpage->text;
    $ctext->font($titletextfont,$titlefontsize/pt);
    $ctext->fillcolor($caxiscolour);
    $ctext->translate(($cplotyaxis+$cplotsize/2)/mm,($cplotxaxis+$cplotsize+5)/mm);
    $ctext->text_center($plottitle." correlation plot");

    $cplot->fillcolor($cbasecolour);
    $cplot->rect($cplotyaxis/mm,$cplotxaxis/mm,$cplotsize/mm,$cplotsize/mm);
    $cplot->fill;
    $cplot->strokecolor($caxiscolour);
    $cplot->rect(($cplotyaxis-1)/mm,($cplotxaxis-1)/mm,($cplotsize+2)/mm,($cplotsize+2)/mm);
    $cplot->stroke;
    for my $y (qw/0 10 20 30 40 50 60 70 80 90 100/){
	my $cypos=$cplotxaxis+$cplotsize*($y-$cplotystart)/$cplotyscale;
	$cplot->move(($cplotyaxis-1)/mm,$cypos/mm);
	$cplot->line(($cplotyaxis-4)/mm,$cypos/mm);
	$cplot->stroke;
	$ctext->font($moltextfont,10/pt);
	$ctext->translate(($cplotyaxis-6)/mm,$cypos/mm - 5/pt);
	$ctext->fillcolor($caxiscolour);
	$ctext->text_right($y);
	#add text labels
    }
    $ctext->font($moltextfont,12/pt);
    $ctext->fillcolor($caxiscolour);
    #print STDERR join(":",$ctext->textpos),"\n";
    $ctext->rotate(90);
    #print STDERR join(":",$ctext->textpos),"\n";
    $ctext->transform(-translate=>[($cplotyaxis-15)/mm,($cplotxaxis+$cplotsize* 0.5 )/mm],-rotate=>90);
    #print STDERR join(":",$ctext->textpos),"\n";
    $ctext->text_center("Genetic Map Position (cM)");

# need ten labels but need to scale appropriately.
    my $maxmolscale=length("".$maxmolpos); #number of digits.
    my $xint="1". "0" x ($maxmolscale-1);
    print STDERR "XINT $xint $maxmolpos\n";
    while ($maxmolpos/$xint <5) {$xint /=2;}
    for (my $xl=0; $xl<$maxmolpos; $xl+=$xint){
	my $xmark=$cplotyaxis+$cplotsize*$xl/$maxmolpos;
	$cplot->move($xmark/mm, ($cplotxaxis-1)/mm);
	$cplot->line($xmark/mm, ($cplotxaxis-4)/mm);
	$cplot->stroke;
	$ctext->font($moltextfont,10/pt);
	$ctext->transform(-translate=>[$xmark/mm - 5/pt, ($cplotxaxis-6)/mm], -rotate=>270);
	$ctext->text($xl/1000);
    }
    $ctext->font($moltextfont,12/pt);
    $ctext->translate(($cplotyaxis+0.5*$cplotsize)/mm, ($cplotxaxis-30)/mm);
    $ctext->text_center("Physical coordinate (kbp)");
    $pdf->finishobjects($cplot,$ctext);
}


# draw chromosome

my $chrob=$page->gfx;
if ($linewidth){
    $chrob->linewidth($linewidth);
}
$chrob->fillcolor($chrcolour);
$chrob->strokecolor($chrcolour);

# draw base chromosome.
#print STDERR join(":",$chrx/mm, $basemargin/mm, ($chrwidth)/mm, ($maxheight)/mm), "\n";
#$chrob->rect($chrx/mm, $basemargin/mm, ($chrwidth)/mm, ($maxheight)/mm);
#$chrob->fill;
#$chrob->circle(($chrx+$chrwidth/2)/mm,$basemargin/mm, ($chrwidth/2)/mm);
#$chrob->fill;
#$chrob->circle(($chrx+$chrwidth/2)/mm,($basemargin+$maxheight)/mm, ($chrwidth/2)/mm);
#$chrob->fill;

# DRAW CHROMOSOME OUTLINE

#draw outline
$chrob->move($chrx/mm, $basemargin/mm);
$chrob->line($chrx/mm,($basemargin+$maxheight)/mm);
$chrob->stroke;
$chrob->move(($chrx+$chrwidth)/mm, $basemargin/mm);
$chrob->line(($chrx+$chrwidth)/mm,($basemargin+$maxheight)/mm);
$chrob->stroke;
#print STDERR join(":",($chrx+$chrwidth/2)/mm,($basemargin+$maxheight)/mm, 270, 90, ($chrwidth/2)/mm,($chrwidth/2)/mm,1 ), "\n";
$chrob->arc(($chrx+$chrwidth/2)/mm,($basemargin+$maxheight)/mm,  ($chrwidth/2)/mm,($chrwidth/2)/mm, 1, 180,1);
$chrob->stroke;
#print STDERR join(":",($chrx+$chrwidth/2)/mm,$basemargin/mm, ($chrwidth/2)/mm,($chrwidth/2)/mm, 181, 360, 1),"\n";
$chrob->arc(($chrx+$chrwidth/2)/mm,$basemargin/mm, ($chrwidth/2)/mm,($chrwidth/2)/mm,181, 360,1);
$chrob->stroke;

#need to draw molecules after doing the layout? or just redraw with scaling?

#DRAW MOLECULES
#print STDERR "plotting molecules\n";
my @molcolfill=qw/white #CCCCCC/;
my $mc=0;
my $molbins=int(($maxheight/mm)/($binsize*$molfontsize/pt));
#calculate label size needed.
if ($molbins < scalar @molecules) {

    $molfontsize=int(($maxheight/mm)/($binsize*(scalar @molecules)/pt));
    $molbins=int(($maxheight/mm)/($binsize*$molfontsize/pt));
#    print STDERR "$molfontsize ",(scalar @molecules)," ",(pt*$maxheight/mm)," ",($molfontsize*$binsize*(scalar @molecules)),"\n";
    if ($molfontsize < $minfont) {
	$molfontsize=$minfont;
	$molbins=int(($maxheight/mm)/($binsize*$molfontsize/pt));
	warn("too many molecules ( ",(scalar @molecules)," in $molbins)");
    }
}
#print STDERR (scalar @molecules)," molecules to plot in $molbins  slots\n";

# now add in binning code for molecule labels.
my @mollabelbins=();

foreach my $m (sort {$a->{start} <=> $b->{start} } @molecules){
    push @{$mollabelbins[int($molbins* ($m->{start}+($m->{length}/2))/$maxmolpos)]}, $m;
    $mc++;
    my $molob=$page->gfx;
    if ($linewidth){
	$molob->linewidth($linewidth);
    }
    $molob->strokecolor($molcolour);
#    print STDERR join(":", $m),"\n";
    my $endy=($m->{length}/$maxmolpos)*$maxheight;
    my $starty=$basemargin+$maxheight-((($m->{start}+$m->{length})/$maxmolpos)* $maxheight); #invert
    if ($invert){
	$starty=$basemargin+(($m->{start}/$maxmolpos)* $maxheight); #invert
    }
    $molob->fillcolor($molcolfill[$mc %2]);
    $molob->rect( $molx/mm, $starty/mm, ($molwidth)/mm, $endy/mm);
    $molob->fill;
    $molob->strokecolor($molcolour);
    $molob->rect( $molx/mm, $starty/mm, ($molwidth)/mm, $endy/mm);
    $molob->stroke;

# draw molecule bands on the correlation and fitted plots. (cpage and dpage)
    my $cmolob=$cpage->gfx;
    $cmolob->fillcolor($molcolfill[$mc %2]);
    $cmolob->rect( ($cplotyaxis+$cplotsize*($m->{start}/$maxmolpos))/mm, $cplotxaxis/mm, $cplotsize * ($m->{length}/$maxmolpos)/mm, $cplotsize/mm);
    $cmolob->fill;
    my $dmolob=$dpage->gfx;
    $dmolob->fillcolor($molcolfill[$mc %2]);
    $dmolob->rect( ($cplotyaxis+$cplotsize*($m->{start}/$maxmolpos))/mm, $cplotxaxis/mm, $cplotsize * ($m->{length}/$maxmolpos)/mm, $cplotsize/mm);
    $dmolob->fill;
    $pdf->finishobjects($molob,$cmolob, $dmolob);
}

#now to draw the labels.
my @labbins=();
my $currentbin=0;
my $lastoccbin=-1;
my $isgroup=0;
my $plotted=0;

# work through the list of labels, plotting in appropriate locations.

while ($currentbin < $molbins){
    
    #print STDERR "current bin $currentbin\n";
    if (!defined( @{$mollabelbins[$currentbin]})){
	$currentbin++; $isgroup=0;next;
#	print STDERR "Bin $currentbin No names - isgroup is $isgroup\n";
    }else{
	#print STDERR "Bin $currentbin occupied - isgroup is $isgroup\n";
	my $plotbin=$lastoccbin+1;
	unless ($isgroup){
	    my $binp=$currentbin;
	    my $binc=0;
	    my $bing=0;
	    while (defined(@{$mollabelbins[$binp]})){
		$bing++;
		$binc+= scalar @{$mollabelbins[$binp]};
	#	print STDERR "Grouping Bin $binp - adding ",(scalar @{$mollabelbins[$binp]})," to give $binc\n";
		$binp++;
	    }
	    # how many locations remain?
	    my $mremain=(scalar @molecules)-$plotted;
	 #   print STDERR "$mremain molecules left to plot in ",($molbins-$lastoccbin)," slots\n";
	  #  print STDERR "at $currentbin. this group of $binc in $bing slots. Last free is $lastoccbin\n"; 
	    while (($molbins-($plotbin) <$mremain) ||
		   (($plotbin>$lastoccbin) &&
		   ($plotbin>($currentbin-int(($binc-$bing)/2))))) {
		$plotbin--;
	    };
	    $plotbin++;
	    $isgroup=1;
	}
	foreach my $ml (@{$mollabelbins[$currentbin]}){
	    plotmol($ml, $plotbin); # plots the text and line for the molecule labels.
	    $lastoccbin=$plotbin;
	    $plotbin++;
	    $plotted++;
	}
    
    }
    $currentbin++;
}

sub plotmol {
    my ($m, $plotbin)=@_;
  #  print STDERR "plotting ",$m->{name}," at $plotbin\n";
 #need to add text here.
    my $moltext=$page->text;
    my $linkob=$page->gfx;
   # print STDERR "plotting with font size $molfontsize\n"; 
    $moltext->font($moltextfont, $molfontsize/pt);
    my $texty=$basemargin+$maxheight - ($maxheight*$plotbin/$molbins);
    my $moly=$basemargin+$maxheight-($maxheight*($m->{start}+($m->{length}/2))/$maxmolpos);
    if ($invert){
	$texty=$basemargin+($maxheight*$plotbin/$molbins);
	$moly=($maxheight*($m->{start}+($m->{length}/2))/$maxmolpos)+$basemargin;
    }
   # print STDERR "text: ",$plotbin/$molbins, "($texty) mol: ", ($m->{start}+($m->{length}/2))/$maxmolpos, "($moly)\n";
    #unless (int($molbins*($m->{start}+($m->{length}/2))/$maxmolpos) == $plotbin){
	#draw line to molecule midpoint
	$linkob->strokecolor(('darkgray','gray')[$plotbin %2]); # molecule label text colour
# alternate between two colours to make it easier to see which molecule is which.
	$linkob->move(($mollabelwidth-2)/mm,($texty/mm)+($molfontsize/2*pt));
	$linkob->line(($mollabelwidth-10)/mm, $moly/mm);
	$linkob->stroke;   
	$linkob->strokecolor('#CCCCCC'); # molecule label connecting line colour.
	$linkob->move(($molx+$molwidth+2)/mm,$moly/mm);
	$linkob->line(($mollabelwidth-10)/mm, $moly/mm);
	$linkob->stroke;   
    
    #}
    $moltext->fillcolor(('darkgray','gray')[$plotbin %2]);
    $moltext->translate(($mollabelwidth)/mm,$texty/mm);
    $moltext->text($m->{name});
    $pdf->finishobjects($linkob, $moltext);

}


#print STDERR "plotted molecules\n";

#DRAW MARKERS

my $markerbins=int(($maxheight/mm)/($binsize*$markerfontsize/pt));
# calculate the number of marker label bins (number of rows) to plot.
if ($molposbin) {
    $markerbins=int($crowd*$maxmolpos/$molposbin)+1;
    $markerfontsize= int(($maxheight/mm)/($binsize*$markerbins/pt));
    if ($markerfontsize <$minfont) {
	$markerfontsize=$minfont;
	$markerbins=int(($maxheight/mm)/($binsize*$markerfontsize/pt));
	$molposbin=int($maxmolpos/($crowd*$markerbins))+1;
	warn("marker plot resolution too high. Bin rescaled to $molposbin\n");
    }
}
print STDERR "Plotting markers into $markerbins bins\n";

# PROCESS MARKERS into bins.

my @seqbins=();

foreach my $mark (@markers) {
    my $b=int($markerbins*$mark->{molpos}/$maxmolpos);
    $mark->{count}=0;

    push @{$seqbins[$b]{$mark->{type}}}, $mark;
    ${$seqbins[$b]{$mark->{type}}}[0]->{count}++;
    #print STDERR "Mark bin $b type $mark->{type} count $mark->{count} bincount ".${$seqbins[$b]{$mark->{type}}}[0]->{count}."\n";
}	
#we can now take just the first marker of each type from the bin.

#foreach my $bin (@seqbins){
#    foreach my $type (sort keys %$bin){
#	for (my $mi=0; $mi<scalar @{$bin->{$type}}; $mi++){
#	    if ($mi){
#		$bin->{$type}[$mi]->{count}=0;
#	    }else{
#		$bin->{$type}[$mi]->{count}=scalar @{$bin->{$type}};
#	    }
#	}
#   }
#}


my @markerbins=();
my @mapbins=();
#print STDERR "markerbins $markerbins\n";
foreach my $m (sort {$a->{molpos} <=> $b->{molpos}?$a->{molpos} <=> $b->{molpos}:$a->{name} cmp $b->{name}} @markers) {
    $m->{bin}=int($markerbins*$m->{molpos}/$maxmolpos);
    if ($m->{count}){
	push @{$markerbins[$m->{bin}]}, $m; 
    } else {
#	print STDERR "Skipping marker $m->{name} as count is zero.\n";
    }
}
foreach my $m (sort {$a->{avemappos} <=> $b->{avemappos}?$a->{avemappos} <=> $b->{avemappos}:$a->{name} cmp $b->{name}} @markers) {
    if ($plotallmap){
	$m->{mapbin}=int($markerbins*$m->{avemappos}/$maxchrpos);
    }else{
	$m->{mapbin}=int($markerbins*binmappos($m->{avemappos})/$maxchrpos);
    }
    #print STDERR "MAP bin ".join(" : ",$m->{mapbin},$markerbins, $m->{avemappos},$maxchrpos, $m->{name}  )."\n";
    if ($m->{mapbin} <0){
	if ($plotallmap){
	    push @{$mapbins[0]{"m_".$m->{avemappos}}}, $m;
	}else{
	    push @{$mapbins[0]{"m_0.0"}}, $m;
	}
    }else{
	if ($plotallmap){
	    push @{$mapbins[$m->{mapbin}]{"m_".$m->{avemappos}}}, $m;
	}else{
	    push @{$mapbins[$m->{mapbin}]{"m_".binmappos($m->{avemappos})}}, $m;
	}
    }
}

foreach my $l (@markers){
#    print STDERR "plotting marker $l->{name}\n";
    my $linkob=$page->gfx;
    my $linktxt=$page->text;
    $linkob->strokecolor($links{$l->{type}}{colour});
    my $clinkob=$cpage->gfx;
    $clinkob->strokecolor($links{$l->{type}}{colour});
## draw cross for marker on correlation plots.
    my $cmy=$cplotxaxis+$cplotsize*($l->{avemappos}-$minchrpos)/$cplotyscale;
    my $cmx=$cplotyaxis+$cplotsize*$l->{molpos}/$maxmolpos;
    $clinkob->move($cmx/mm, ($cmy-2)/mm);
    $clinkob->line($cmx/mm, ($cmy+2)/mm);
    $clinkob->stroke;
    $clinkob->move(($cmx-2)/mm, $cmy/mm);
    $clinkob->line(($cmx+2)/mm,$cmy/mm);
    $clinkob->stroke;
    my $dlinkob=$dpage->gfx;
    $dlinkob->strokecolor($links{$l->{type}}{colour});
    $dlinkob->move($cmx/mm, ($cmy-2)/mm);
    $dlinkob->line($cmx/mm, ($cmy+2)/mm);
    $dlinkob->stroke;
    $dlinkob->move(($cmx-2)/mm, $cmy/mm);
    $dlinkob->line(($cmx+2)/mm,$cmy/mm);
    $dlinkob->stroke;
# need to check for a range of positions
    my $lcy=0;
    my $lmy=$basemargin+$maxheight-($maxheight*$l->{molpos}/$maxmolpos);
    if ($invert){
	$lmy=($maxheight*$l->{molpos}/$maxmolpos)+$basemargin;
    }
# draw marker line
    if ($l->{mappos} =~ m/(\d+)-(\d+)/) {
	my $s=$1;
	my $e=$2;
	my $m=$l->{avemappos};
	my $lcs=$basemargin+$maxheight-($maxheight*$s/$maxchrpos);
	if ($invert) {
	    $lcs=($maxheight*$s/$maxchrpos)+$basemargin;
	}
	$lcy=$basemargin+$maxheight-($maxheight*$m/$maxchrpos);
	if ($invert){
	    $lcy=($maxheight*$m/$maxchrpos)+$basemargin;
	}
	$linkob->move(($chrx)/mm,$lcs/mm);
	$linkob->line(($chrx+$chrwidth)/mm, $lcs/mm);
	$linkob->stroke;
	my $lce=$basemargin+$maxheight-($maxheight*$e/$maxchrpos);
	if ($invert){
	    $lce=($maxheight*$e/$maxchrpos)+$basemargin;
	}
	$linkob->move(($chrx)/mm,$lce/mm);
	$linkob->line(($chrx+$chrwidth)/mm,$lce/mm);
	$linkob->stroke;
	$linkob->move(($chrx+$chrwidth)/mm, $lcs/mm);
	$linkob->line(($chrx+$chrwidth+2)/mm, $lcy/mm);
	$linkob->line(($chrx+$chrwidth)/mm, $lce/mm);
	$linkob->stroke;
	$linkob->move(($chrx)/mm, $lcs/mm);
	$linkob->line(($chrx-2)/mm, $lcy/mm);
	$linkob->line(($chrx)/mm, $lce/mm);
	$linkob->stroke;

    }else {
	$lcy=$maxheight+$basemargin-($maxheight*$l->{mappos}/$maxchrpos);
	if ($invert){
	    $lcy=($maxheight*$l->{mappos}/$maxchrpos)+$basemargin;
	}
	$linkob->move(($chrx-2)/mm,$lcy/mm);
	$linkob->hline(($chrx+$chrwidth+2)/mm);
	$linkob->stroke;
    }
    $linkob->move(($molx-2)/mm,$lmy/mm);
    $linkob->hline(($molx+$molwidth+2)/mm);
    $linkob->stroke;

    #if ($links{$l->{type}}{line} eq 'dash'){
#	$linkob->linedash(2,2);
    #} else {
#	$linkob->linedash(0);
 #   } 
# draw molecule side horiz line
    $linkob->move(($molx-2)/mm,$lmy/mm);
    $linkob->line(($chrx+2+$chrwidth)/mm,$lcy/mm);
    $linkob->stroke;
 #   print STDERR "font $links{$l->{type}}{font}\n";
    $pdf->finishobjects($linkob,$clinkob,$dlinkob);
}

#$pdf->end();


#label layout.

#Start at top occupied bin.
#1 entity -> plot.
#Is there sufficient space above or below to plot? Plot.
#Plot first n and last m in the n and m spaces above and below. Plot remainder spread in row 2 filling up first.


$currentbin=0;
my $lastoccbin1=-1;
my $lastoccbin2=-1;

while ($currentbin < $markerbins){
    
    # TODO - plot number of markers in a bin and the bin value, not individuals.

    #print STDERR "current bin $currentbin\n";
    if (!defined( @{$markerbins[$currentbin]})){
	$currentbin++; next;
    } elsif (scalar @{$markerbins[$currentbin]}==1 ){
	my $plotbin=$currentbin;
	while ($plotbin<=$lastoccbin1) {$plotbin++;}
	plotmarker($markerbins[$currentbin]->[0],$plotbin,0);
	$lastoccbin1=$plotbin;
	$currentbin++;
    }else {
	my $toplot=scalar @{$markerbins[$currentbin]};
	my $plotbin=$currentbin;
	while ($plotbin<=$lastoccbin1) {$plotbin++;}
	my $nextoccbin=$plotbin+1;
	while (!defined(@{$markerbins[$nextoccbin]}) && $nextoccbin < $markerbins){
	    $nextoccbin++;
	}
	

	if ($nextoccbin-$lastoccbin1 >$toplot) {
	    for (my $b=0; $b<$toplot; $b++){
		my $plotbin=1+$lastoccbin1;
		while ($plotbin < $currentbin-$toplot){ $plotbin++;}
		plotmarker($markerbins[$currentbin]->[$b], $plotbin,0);
		$lastoccbin1=$plotbin;
	    }
	}else{
	    my $upset=$currentbin-($lastoccbin1+1);
	    my $upset2=$toplot-$upset;
	    my $start2=$currentbin-$upset2;
	    while ($start2<=$lastoccbin2) {
		$start2++;
	    }
	    for (my $p=0; $p<$upset; $p++){
		my $plotbin=$currentbin+$p-$upset;
		plotmarker($markerbins[$currentbin]->[$p],$plotbin,0);
		$lastoccbin1=$plotbin;
	    }
	    $lastoccbin1=$currentbin;
	    for (my $s=0; $s<$toplot-$upset;$s++){
		my $plotbin=$start2+$s;
		plotmarker($markerbins[$currentbin]->[$s+$upset],$plotbin, 1);
		$lastoccbin2=$plotbin;
	    }
	}
	$currentbin++;
    }
}


## Plotting genetic map markers
# plot unanchored
foreach my $m (@unanchored){plotunmap($m);}
#plot anchored
$currentbin=0;
$lastoccbin1=-1;
$lastoccbin2=-1;
while ($currentbin <= $markerbins){
    my @keys=();
    #print STDERR "current bin $currentbin\n";
    if (!defined( %{$mapbins[$currentbin]})){
	$currentbin++; next;
    } else {
	@keys=sort {substr($a,2) <=> substr($b,2)} keys %{$mapbins[$currentbin]};
    }
    if (scalar @keys==1 ){
	my $plotbin=$currentbin;
	while ($plotbin<=$lastoccbin1) {$plotbin++;}
	
	plotmap($keys[0], scalar @{$mapbins[$currentbin]{$keys[0]}}, $plotbin,0);

	$lastoccbin1=$plotbin;
	$currentbin++;
    }else {
	my $toplot=scalar @keys;
	my $plotbin=$currentbin;
	while ($plotbin<=$lastoccbin1) {$plotbin++;}
	my $nextoccbin=$plotbin+1;
	while (!defined(%{$mapbins[$nextoccbin]}) && $nextoccbin < $markerbins){
	    $nextoccbin++;
	}
	if ($nextoccbin-$lastoccbin1 >$toplot) {
	    for (my $b=0; $b<$toplot; $b++){
		my $plotbin=1+$lastoccbin1;
		while ($plotbin < $currentbin-$toplot){ $plotbin++;}
		plotmap($keys[$b], scalar @{$mapbins[$currentbin]{$keys[$b]}}, $plotbin,0);
		$lastoccbin1=$plotbin;
	    }
	}else{
	    my $upset=$currentbin-($lastoccbin1+1);
	    my $upset2=$toplot-$upset;
		my $start2=$currentbin-$upset2;
	    while ($start2<=$lastoccbin2) {
		$start2++;
	    }
	    for (my $p=0; $p<$upset; $p++){
		my $plotbin=$currentbin+$p-$upset;
		plotmap($keys[$p],scalar @{$mapbins[$currentbin]{$keys[$p]}},$plotbin,0);
		$lastoccbin1=$plotbin;
	    }
	    $lastoccbin1=$currentbin;
	    for (my $s=0; $s<$toplot-$upset;$s++){
		my $plotbin=$start2+$s;
		plotmap($keys[$s+$upset], scalar @{$mapbins[$currentbin]{$keys[$s+$upset]}},$plotbin, 1);
		    $lastoccbin2=$plotbin;
	    }
	}
	$currentbin++;
    }
}


my $titletxt=$page->text;
$titletxt->font($titletextfont,$titlefontsize/pt);
$titletxt->translate(($molx+$molwidth/2)/mm, ($basemargin+$maxheight+5)/mm);
$titletxt->fillcolor($titlecolour);
my $fulltitle=$plottitle;
unless ($nosize){
    $fulltitle .= " ($maxmolpos bp)";
}
$titletxt->text_center($fulltitle);
$pdf->finishobjects($titletxt);

# need to pull back all markers in rearranged molecule and plot
my @repmarkers=sort {int($a->{molpos} <=>$b->{molpos})} $assembly->getMarkers($mapname);
print STDERR "Retrieved ".(scalar @repmarkers)." markers\n";
foreach my $m (@repmarkers){
    my $elinkob=$dpage->gfx;
    $elinkob->strokecolor($links{$m->{type}}{colour});
    my $cmy=$cplotxaxis+$cplotsize*($m->{mappos}-$minchrpos)/$cplotyscale;
    my $cmx=$cplotyaxis+$cplotsize*$m->{molpos}/$maxmolpos;
    
    $elinkob->circle($cmx/mm, $cmy/mm, 1/mm);
    $elinkob->stroke;
    $pdf->finishobjects($elinkob);
}
my $oldline=$dpage->gfx;
$oldline->strokecolor('red');
my $sp=$repmarkers[0]->{molpos};
my $ep=$repmarkers[$#repmarkers]->{molpos};
#print STDERR "plotting original curve between $sp and $ep\n";
my $step=int(($ep-$sp)/200);# ca. 200 steps to draw the line.
my $lastx=0;
my $lastyold=0;
#print "parameters are $report->{assembly}{final}{a0} $report->{assembly}{final}{a1} $report->{assembly}{final}{a2} $report->{assembly}{final}{a3}\n";
for (my $pos=$sp; $pos<$ep; $pos += $step){
    my $xc=$cplotyaxis+$cplotsize*$pos/$maxmolpos;
    my $yo=gety($pos, $report,0);
    my $yco=$cplotxaxis+$cplotsize*($yo-$minchrpos)/$cplotyscale;
    #print "plotting coordinates $xc, $yco\n";
    if ($lastx){
	$oldline->move($lastx/mm,$lastyold/mm);
	$oldline->line($xc/mm,$yco/mm);
	$oldline->stroke;
    }else{
	$oldline->move($xc/mm,$yco/mm);
    }
    $lastx=$xc;
    $lastyold=$yco;
}
my $newline=$dpage->gfx;
$newline->strokecolor('blue');
 $lastx=0;
my $lastynew=0;
for (my $pos=$sp; $pos<$ep; $pos += $step){
    my $xc=$cplotyaxis+$cplotsize*$pos/$maxmolpos;
    my $yn=gety($pos, $report,1);
    my $ycn=$cplotxaxis+$cplotsize*($yn-$minchrpos)/$cplotyscale;
    if ($lastx){
	$newline->move($lastx/mm, $lastynew/mm);
	$newline->line($xc/mm,$ycn/mm);
	$newline->stroke;
    }else{
	$newline->move($xc/mm,$ycn/mm);
    }
    $lastx=$xc;
    $lastynew=$ycn;
}
$pdf->finishobjects($oldline,$newline);
#now plot the report as a table.




my $pdftable=PDF::Table->new();

my %headerprops =(
    font       => $font{Helvetica}{Bold},
    font_size  => 9,
    font_color => '#006666',
    bg_color   => 'lightyellow',
    repeat     => 1  
);
my @cellprops=(); # array to capture detailed cell formatting.
my @colprops=(); # array of hashrefs for column formatting
my @tabledata=(); #plot data = array of arrays

my $margin=15; # page margin. Page is $width*$height

my $rtext=$rpage->text();
$rtext->font($titletextfont,$titlefontsize/pt);
my $curry=($height-$margin)/mm-(1.3*$titlefontsize/pt);
$rtext->translate(($margin)/mm,$curry);
$rtext->fillcolor($titlecolour);
$rtext->text("Report for $plottitle");
$pdf->finishobjects($rtext);

$curry -= 5/mm;
#add headers to table
print STDERR "FIT table\n";
push @tabledata,['', "Mean Square Error", "Mean Error", qw/a b c d/];
#print STDERR "$report\n";
push @tabledata,["Initial Fit", sprintf("%5.3e", $report->{assembly}{original}{chi2} ), sprintf("%5.3e", $report->{assembly}{original}{var}), sprintf("%5.3e", $report->{assembly}{original}{a0}), sprintf("%5.3e", $report->{assembly}{original}{a1}), sprintf("%5.3e", $report->{assembly}{original}{a2}), sprintf("%5.3e", $report->{assembly}{original}{a3})];

push @tabledata,["Optimal Fit", sprintf("%5.3e", $report->{assembly}{final}{chi2}), sprintf("%5.3e", $report->{assembly}{final}{var}),sprintf("%5.3e",  $report->{assembly}{final}{a0}), sprintf("%5.3e", $report->{assembly}{final}{a1}), sprintf("%5.3e", $report->{assembly}{final}{a2}), sprintf("%5.3e", $report->{assembly}{final}{a3})];

# need to draw the fit lines.



my $pagecount=0;
#print STDERR "Header table ", $curry * mm, "\n";
($rpage,$pagecount, $curry)=$pdftable->table($pdf, $rpage,\@tabledata, 
					     x=>$margin/mm,
					     w=>($width-(2*$margin))/mm, 
					     start_y=> $curry, 
					     start_h=>100/mm,
					     padding=>5, 
					     font=> $font{Helvetica}{Bold},
					     font_size      => 9,
					     font_color_odd => 'blue',
					     font_color_even=> 'black',
					     background_color_odd  => 'lightgray',         #cell background color for odd rows
					     background_color_even => 'lightblue',     #cell background color for even rows
					     header_props   => \%headerprops
    );

#print STDERR "Header table ", $curry * mm, "\n";
$curry -=5/mm;
if ($curry <75/mm){
    $rpage=$pdf->page;
    $curry=($height-$margin)/mm;
}
 # %{$report{assembly}{original}}=$self->fit();
 #   %{$report{assembly}{final}}=$self->bestfit();
print STDERR "Scaffold table\n";
@tabledata=(["Scaffold", "Start", "Length", "Markers", "Strand", "Mean Square Error", "Mean Error"] );
@cellprops=([{},{},{},{},{},{},{}]);
foreach my $m (sort {$a->{start} <=> $b->{start}} @{$report->{molecules}}){
    my @cp=({},{},{},{},{},{},{});
    my @cn=({},{},{},{},{},{},{});
    my @td=($m->{name}, $m->{start}, $m->{length},$m->{markers});
    my @nd=();
    if ($m->{markers}==0){
	push @td, qw/- - -/;
	for (my $i=0; $i< scalar @cp; $i++){
	    $cp[$i]={background_color=>'lightgray'};
	}
    } elsif ($m->{markers}==1){
	push @td, $m->{originalorient}, sprintf("%5.3e", $m->{originalchi2}), sprintf("%5.3e", $m->{originalvar});
    } else {
	push @td, $m->{originalorient}, sprintf("%5.3e", $m->{originalchi2}), sprintf("%5.3e", $m->{originalvar});
	if ($m->{originalorient} ne $m->{neworient}){
	    @nd=("", "", "inverted","", $m->{neworient}, sprintf("%5.3e", $m->{newchi2}), sprintf("%5.3e", $m->{newvar}));
	    $cn[4]{font_color}='red';
	    $cp[0]{background_color}='Tomato';
	} else {
	    $cp[0]{background_color}='PaleGreen';
	}
	
	if ($m->{originalchi2} < $m->{newchi2}){
	    $cp[5]{font_color}='red';
	} else {
	    $cn[5]{font_color}='red';
	}
	if ($m->{originalvar} < $m->{newvar}){
	    $cp[6]{font_color}='red';
	} else {
	    $cn[9]{font_color}='red';
	}
    }
   # print STDERR join(":","Moldata", @td),"\n";
    push @tabledata, \@td;
    push @cellprops, \@cp;
    if (scalar @nd >2){ push @tabledata, \@nd;    push @cellprops, \@cn; }
}


 #   @{$report{molecules}}=();
 #   foreach my $mol (sort {$self->{molecules}{$a}{start} <=> $self->{molecules}{$b}{start}} keys %{$self->{molecules}}){
#	my %m=%{$self->{molecules}{$mol}};
#	my $mstats={
#	    name=>$m{name},
#	    start=>$m{start},
#	    originalorient=>$m{originalorient}?"reverse":"forward",
#	};
#	if (exists($m{chi2})){
#	    $mstats->{originalchi2}=$m{chi2};
#	    $mstats->{originalvar}=$m{var};
#	    $mstats->{neworient}=$m{molecule}->isFlipped()?"reverse":"forward";
#	    my %newfit=$self->_fitMolecule($mol);
##	    $mstats->{newchi2}=$newfit{chi2};
#	    $mstats->{newvar}=$newfit{var};
#	    $mstats->{markers}=$newfit{markers};

my @shortdata=@tabledata[0..1];
my @shortcell=@cellprops[0..1];

#print STDERR (scalar @tabledata),"   ",(scalar @cellprops),"\n";

($rpage,$pagecount, $curry)=$pdftable->table($pdf, $rpage,\@tabledata, 
					     x=>$margin/mm,
					     w=>($width-(2*$margin))/mm, 
					     start_y=> $curry, 
					     next_y=>($height-2*$margin)/mm,
					     start_h=> $curry-$margin/mm, 
					     next_h=>($height-3*$margin)/mm,
					     padding=>5, 
					     font=> $font{Helvetics}{Bold},
					     font_size      => 8,
					     font_color=> 'black',
					     header_props   => \%headerprops,
					     #cell_props =>\@shortcell,
    );

#print STDERR "$rpage, $pagecount, ",($curry * mm),"\n";
$curry=$curry-5/mm;

# calculate mean marker map position and rank molecules by that.

my %molrank =();
my %avemap=();
my $rank=0;
foreach my $m (sort {$a->{start} <=>$b->{start}} @molecules){
	if (0){
    	warn "Molecule listed:\n";
    	foreach  my $k (keys %$m){
    		warn "$k: $m->{$k}\n";
    	}
    			
    }
    
	
    my %mol=$assembly->getMolecule($m->{name}, $m->{fragment});
    if (0){
    	warn "Molecule returned:\n";
    	foreach  my $k (keys %mol){
    		warn "$k: $mol{$k}\n";
    	}
    			
    }
    if ((scalar $mol{molecule}->getMarkers($mapname))>0){
	$rank++;
	$molrank{join("_",$m->{name}, $m->{fragment})}=$rank;
	$avemap{join("_",$m->{name}, $m->{fragment})}=$mol{molecule}->aveMapPos($mapname);
    }
}
my %maprank=();
$rank=0;
foreach my $m (sort { $avemap{$a} <=> $avemap{$b} } keys %molrank){
    $rank++;
    $maprank{$m}=$rank;
}

print STDERR "RANK table\n";

my @ranktable=(["Molecule", "Assembly order", "Map order", "Difference"]);

foreach my $m (sort { $molrank{$a} <=> $molrank{$b} } keys %molrank){
    my $diff=$maprank{$m}-$molrank{$m};
    my $difftext="$diff";
    if ($diff==0){ $difftext="-";}

    push @ranktable, [$m,$molrank{$m}, $maprank{$m}, $difftext];
}
if ($debug){
	foreach my $r (@ranktable){
		print STDERR join(":","MAPRANK", @$r)."\n";
	}
}

if ($curry <75/mm){
    $rpage=$pdf->page;
    $curry=($height-$margin)/mm;
}

($rpage,$pagecount, $curry)=$pdftable->table($pdf, $rpage,\@ranktable, 
					     x=>$margin/mm,
					     w=>($width-(2*$margin))/mm, 
					     start_y=> $curry-($margin/mm), 
					     next_y=>($height-2*$margin)/mm,
					     start_h=> $curry-$margin/mm, 
					     next_h=>($height-3*$margin)/mm,
					     padding=>5, 
					     font=> $font{Helvetica}{Bold},
					     font_size      => 10,
					     font_color=> 'black',
					     header_props   => \%headerprops,
    );

$curry -= 5/mm;

if ($curry < 75/mm) {
    $rpage=$pdf->page;
    $curry=($height-2*$margin)/mm;
}

# pick up markers and order by error.

my %markererror=();

foreach my $m (@markers){ 

    my $y =gety($m->{molpos},$report,0);
    $markererror{$m->{name}}={error=>abs($y-$m->{avemappos}), molpos=>$m->{molpos}, calcmappos=>$y, givenmappos=>$m->{avemappos}};
}

print STDERR "MARKER table\n";
my @errortable=(["Marker", "sequence position", "map position", "calculated fit", "error"]);

foreach my $m (sort {int($markererror{$b}{error} <=>$markererror{$a}{error})} keys %markererror){
    push @errortable, [ $m, $markererror{$m}{molpos}."bp", $markererror{$m}{givenmappos}, sprintf("%5.3f",$markererror{$m}{calcmappos}),sprintf("%5.3f", $markererror{$m}{error})];
#    print STDERR join(":",$m, $markererror{$m}{molpos}."bp", $markererror{$m}{givenmappos}, sprintf("%5.3f", $markererror{$m}{calcmappos}), sprintf("%5.3f",$markererror{$m}{error})),"\n";
}
#print STDERR (scalar @errortable),"\n"; 

my @plotlist=();
my $plotlen=50;
if (scalar @errortable <$plotlen){
    $plotlen=scalar @errortable-1;
}
foreach my $e (@errortable[0..$plotlen]){
    if ($e){
#	print STDERR join("::", @$e),"\n";
	push @plotlist, $e;
    } 

}
if ($plotlist[0]){

    ($rpage,$pagecount, $curry)=$pdftable->table($pdf, $rpage,\@plotlist,
					     x=>$margin/mm,
					     w=>($width-(2*$margin))/mm, 
					     start_y=> $curry, 
					     next_y=>($height-2*$margin)/mm,
					     start_h=> $curry-$margin/mm, 
					     next_h=>($height-2*$margin)/mm,
					     padding=>5, 
					     font=> $font{Helvetica}{Bold},
					     font_size      => 10,
					     font_color=> 'black',
					     header_props   => \%headerprops,
    );
}




sub plotmarker {
    my ($l, $plotbin,$column)=@_;
 #   print STDERR "plotting marker $l in $plotbin column $column\n";
    my $linkob=$page->gfx;
    $linkob->strokecolor($links{$l->{type}}{colour});
    my $linktxt=$page->text;
    my ($ff, $fs)=split / /, $links{$l->{type}}{font};
  #  print STDERR "font $font{$ff}{$fs}\n";
    my $labelpos=$basemargin+$maxheight-($maxheight*($plotbin+0.5)/$markerbins);
    my $lcy=$basemargin+$maxheight-($maxheight*$l->{avemappos}/$maxchrpos);
    my $lmy=$basemargin+$maxheight-($maxheight*$l->{molpos}/$maxmolpos);
    my $sb=int($l->{molpos}/$maxmolpos);
    my $mc=exists($seqbins[$sb]{$l->{type}})?scalar (@{$seqbins[$sb]{$l->{type}}}):0;
    if ($invert){
	$labelpos=($maxheight*($plotbin+0.5)/$markerbins)+$basemargin;
	$lcy=($maxheight*$l->{avemappos}/$maxchrpos)+$basemargin;
	$lmy=($maxheight*$l->{molpos}/$maxmolpos)+$basemargin;
    }
    if ($l->{mappos} =~ m/(\d+)-(\d+)/) {
	my $s=$1;
	my $e=$2;
	my $m=$l->{avemappos};
	$lcy=($maxheight*$m/$maxchrpos)+$basemargin;
	    }
    if (exists($font{$ff}{$fs})) {
	$linktxt->font($font{$ff}{$fs}, $markerfontsize/pt);
    }else {
	$linktxt->font($font{Helvetica}{Roman}, $markerfontsize/pt);
    }
    $linktxt->fillcolor($links{$l->{type}}{colour});


# Plot marker name 

   my $labelx=($molx+$molwidth+2);
    $linkob->move($labelx/mm,$lmy/mm);
    if ($column){
	$labelx+=$markerlabelwidth;
	$linkob->line($labelx/mm,$lmy/mm);
	$linkob->stroke;
	$linkob->move($labelx/mm,$lmy/mm);
    }
    $linkob->line(($labelx+$lineanglespace)/mm,$labelpos/mm);
    $linkob->stroke;
    
    $linktxt->translate(($labelx+$lineanglespace+2)/mm, (($labelpos)/mm) -(($markerfontsize/2)/pt));
    my $labeltext=$l->{name};
    if ($mc >1) {
	$labeltext .= " ($mc)";
    }
    $linktxt->text($labeltext);

    $pdf->finishobjects($linkob, $linktxt);
}

sub plotunmap {
    my %marker=%{$_[0]};
    my $ypos=$basemargin +$maxheight -($maxheight*$marker{'pos'}/$maxchrpos);
    if ($invert) {
	$ypos=$basemargin + ($maxheight*$marker{'pos'}/$maxchrpos);
    }
    my $go=$page->gfx;
    $go->strokecolor('gray');
    my $labelx=($chrx);
    $go->move($labelx/mm,$ypos/mm);
    $labelx+=$chrwidth;
    $go->line($labelx/mm,$ypos/mm);
    $go->stroke;	
    $pdf->finishobjects($go);
}

 sub plotmap {
    my ($l,$count, $plotbin,$column)=@_;
    my $lval=$l;
    $lval=~s/m_//;
    #print STDERR "plotting marker $l in $plotbin column $column\n";
  my $linkob=$page->gfx;
    $linkob->strokecolor('black');
    my $linktxt=$page->text;
  #  print STDERR "font $font{$ff}{$fs}\n";
    my $labelpos=$basemargin+$maxheight-($maxheight*($plotbin+0.5)/$markerbins);
    if ($invert){
	$labelpos=($maxheight*($plotbin+0.5)/$markerbins)+$basemargin;
    }
    my $lcy="";
    if ($lval=~ m/(\d+)-(\d+)/) {
	my $s=$1;
	my $e=$2;
	my $m=($s+$e)/2;
	$lcy=$basemargin+$maxheight-($maxheight*$m/$maxchrpos);
	if ($invert){
	    $lcy=($maxheight*$m/$maxchrpos)+$basemargin;
	}
    } else {
	$lcy=$basemargin+$maxheight-($maxheight*$lval/$maxchrpos);
	if ($invert){ 
	    $lcy=($maxheight*$lval/$maxchrpos)+$basemargin;
	}
    }
    $linktxt->font($font{Helvetica}{Roman}, $markerfontsize/pt);
    my $linkfillcol='black';
    if ($count>1) { $linkfillcol='red';}
    $linktxt->fillcolor($linkfillcol);


# Plot marker name 

   my $labelx=($chrx -2);
    $linkob->move($labelx/mm,$lcy/mm);
    if ($column){
	$labelx-=$maplabelwidth;
	$linkob->line($labelx/mm,$lcy/mm);
	$linkob->stroke;
	$linkob->move($labelx/mm,$lcy/mm);
    }
    $linkob->line(($labelx-$lineanglespace)/mm,$labelpos/mm);
    $linkob->stroke;
    
    $linktxt->translate(($labelx-($lineanglespace+2))/mm, (($labelpos)/mm) -(($markerfontsize/2)/pt));
    my $plural=$count>1?"s":"";
    my $marktext="[$count] $lval";
    if ($nocount){
	$marktext="$lval";
    }
    $linktxt->text_right($marktext);

    $pdf->finishobjects($linkob, $linktxt);

}

$pdf->save;#as($filename);

sub gety {
    my ($pos, $report, $new)=@_;
    my $fit=$new?"final":"original";
    #my $y=$report->{assembly}{$fit}{a0} + $pos * $report->{assembly}{$fit}{a1} + $pos * $pos * $report->{assembly}{$fit}{a2} + $pos *$pos * $pos * $report->{assembly}{$fit}{a3} ;
    my $y=$report->{assembly}{$fit}{a0} + $report->{assembly}{$fit}{a1} * sinh( ($pos - $report->{assembly}{$fit}{a2}*$report->{assembly}{$fit}{size})*$report->{assembly}{$fit}{a3}/$report->{assembly}{$fit}{size});
#    print STDERR "getting calculated molpos for $pos : $y\n";
    return $y;
}

sub binmappos {
    my $pos=shift;
    my $newpos=0;
    if ($pos=~m/\d-\d/){
	my @mp=split/-/, $pos,2;
	my @nmp=();
	foreach my $p (@mp){
	    push @nmp, sprintf("%.${mapdp}f", int($p/$genbin)*$genbin);
	}
	$newpos=join("-",@nmp);
    } else {
	$newpos=sprintf("%.${mapdp}f", int($pos/$genbin)*$genbin);
    }
    return $newpos;
}
sub avemappos {
    my $pos=shift;
    my ($f, $l)=split/-/, $pos,2;
    my $p=$f;
    unless($p){
	$p=-$l;
    }
    if ($l && $f){
	$p=($f+$l)/2;
    }
    my $newpos=sprintf("%.${mapdp}f", int($p/$genbin)*$genbin);
    return $newpos;
}
