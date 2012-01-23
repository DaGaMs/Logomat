# This is HMM::Profile written by Benjamin Schuster-Boeckler

package HMM::Profile;

=head1 NAME

HMM::Profile - Class representing Profile-HMMs

=head1 SYNOPSIS

use HMM::Profile;


my $file = $ARGV[0];    #read hmm-filename from command line

my $pHMM = HMM::Profile->new(-hmmerfile=>$file) || die("Couldn't open $file!\n"); #try to create a new hmm-object from the file

$pHMM->draw_logo(-file      => 'outfile.png',
                -xsize      => $pHMM->length() * 30,
                -ysize      => 400,
                -x_title    => 'HMM-States',
                -y_title    => 'Inf-Content',
                -greyscale  => 1,
                -title      => $pHMM->name(),
                -height_logodds => 0
                );

=head1 DESCRIPTION

HMM::Profile has the ability to draw a Logo from a HMMer file by calling the draw_logo method.

=head1 METHODS

=cut

use strict;
use warnings;
use vars '@ISA'; #inheritance
use LWP::UserAgent;
use HMM;
use HMM::Utilities;
use HMM::Utilities3;

@ISA = qw(HMM); #inherit from HMM

our $REGFONT = '/Library/Fonts/Arial Unicode.ttf';
our $BOLDFONT = '/Library/Fonts/Arial Unicode.ttf';

use PDL::LiteF;
use Imager ':handy';
use Imager::Fill;


=head2 new(%args)

 Usage      : $pHMM = HMM::Profile->new(%args)
 Function   : construct a new HMM object
 Returns    : true, if creation succeeded, else false
 Args       :
    -hmmerfile      # A filename or -handle, HMMer plan7 format
    -pfamid         # The name of a pfam profile. Will be retrieved automatically through LWP
    # The following options will only be parsed if no file or pfamid is given
    -alphabet       # An array reference of letters
    -emissions      # A 3-dimensional PDL matrix, containing both Insert and Match emission probabilities
    -transitions    # A 2-dimensional PDL matrix of transition probabilities
    -nullEmissions  # background probabilities; defaults to uniform distribution
    -startTransitions   # a PDL vector with 2 values; defaults to pdl(1, 0)
    -accession      # An accession for the entry; optional
    -version        # Version; optional
    -name           # Name; optional
    -description    # Description; optional
=cut

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless ($self, ref($class) || $class);

    my $file;

    if($args{-hmmerfile} and ref($args{-hmmerfile}) eq 'Fh' or ref($args{-hmmerfile}) eq 'GLOB' or ref($args{-hmmerfile}) eq 'IO::File'){
        my $fh = $args{-hmmerfile};
        undef $/;
        $file = <$fh>;
        $/ = "\n";
    }
    elsif($args{-hmmerfile}){
        open(HMM, $args{-hmmerfile}) || warn "Couldn't open ".$args{-hmmerfile}.": $!\n" && return 0;
        undef $/;
        $file = <HMM>;
        $/ = "\n";
        close(HMM);
    }
    elsif($args{-pfamid}){
        use LWP::UserAgent;

        # Create a user agent object
        my $ua = LWP::UserAgent->new();
        $ua->agent('HMM.pm/'.$HMM::VERSION);

        # Pass request to the user agent and get a response back
        my $res = $ua->get('http://pfam.sanger.ac.uk/family/hmm?entry='.$args{-pfamid}.'&mode=ls');

        # Check the outcome of the response
        if ($res->is_success) {
            $file = $res->content();
        }
        else { warn "LWP failed! Couldn't get HMM file from Pfam!\n"; return 0; }
    }
    else
    {
        warn("No HMM file or Pfam id passed. Initialising as an HMM object");
        $self = HMM->new(%args);
        bless $self, ref $class || $class;
    }
    if ($file)
    {
        #remove xhtml - encapsulation sent with pfam-profiles
        $file =~ s/<\/?pre>//gi;
        #remove windows line breaks
        $file =~ s/\r\n/\n/g;
        # Parse the *.hmm record
        if ($file =~/^HMMER(\S+)/) {
            $self->{version} = $1;
            if($self->version =~ /^2/){
                $self->_parseFile($file);
            }
            elsif ($self->version =~ /^3/)
            {
                $self->_parseFile3($file);
            }
            else
            {
                warn("Unknown HMMER version used: ", $self->version);
                return;
            }
        }
        else
        {
            warn("First line does not contain HMMER version description, aborting");
            return;
        }
    }
    return $self;
}

sub _raw_data {
    my $self = shift;
    my $height_logodds = shift;
    #initialize the most important variables first

    my $ICM = toICM($self->emissions, $self->nullEmissions, $height_logodds);
    my $maxInfContent = max($ICM->xchg(0,1)->sumover());

    #Width
    my $HPM = toHPM($self->startTransitions, $self->transitions); # calculate Hitting Prob.
    my $motifSize = $ICM->getdim(0);
    # scale the Insert states relative to the estimated retention time, which is 1/(1-p(I->I)) = 1/p(I->M)
    my $width = zeroes($motifSize);
    $width->slice("0:".($motifSize-2).":2") .= $HPM->slice(':,0');
    $width->slice("1:".($motifSize-1).":2") .= ($HPM->slice(':,1')/$self->{'transitions'}->slice(':,3'))->badmask(0); #scale Insert-width by estimated retention period



    return {"width" => $width, "ICM" => $ICM, "maxIC" => $maxInfContent, "HPM" => $HPM};
}

sub flatten {
    my $self = shift;
    my $height_logodds = shift || 0;
    my $data = $self->_raw_data($height_logodds);


    my $ICM = $data->{ICM};
    my $maxInfContent = $data->{maxIC};
    my $HPM = $data->{HPM};
    my $width = $data->{width};

    my @tmp;
    # convert to flat arrays
    for (my $i=0; $i<$ICM->getdim(0); $i++) {
        $tmp[$i] = [$ICM->slice("$i,:")->list];
    }
    $ICM = [@tmp];
    @tmp = ();
    # convert to flat arrays
    for (my $i=0; $i<$HPM->getdim(0); $i++) {
        $tmp[$i] = [$HPM->slice("$i,:")->list];
    }
    $HPM = \@tmp;

    return {"width" =>  [ $width->list ],
            "ICM" => $ICM,
            "maxIC" => $maxInfContent,
            "HPM" => $HPM
    };
}

=head2 draw_logo()

 Usage      : $pHMM->draw_logo(%args)
 Function   : draw a logo from the represented HMM
 Returns    : true, if Image written properly, else false
 Args       :

    -file,           # name of the output file; should end in .png
    -graph_title,    # title; none if empty
    -x_title,        # x-axis descritption; none if empty
    -y_title,        # y-axis description; none if empty
    -xsize,          # Image width; default is 600
    -ysize           # Image heigth; default is 360
    -startpos        # offset-position in the profile
    -endpos          # end-position in the profile
    -greyscale       # create a colorless picture
    -height_emission # divide information content height according to emission probabilities (default);
    -height_logodds  # by default, information content height is divided according to emission probabilities;
                     # if non-zero, heights will be based on log-odds scores (only positive scoring residues)
    -regular_font    # provide full path to true type font file
    -bold_font       # provide full path to true type font file for bold text
    -letter_colours  # reference to hash that asssigns hex colour to residue char
    -title_font_size
    -axis_font_size

=cut

sub draw_logo {
    my $self = shift;

    #default arguments
    my %args = (-xsize      => 600,
                -ysize        => 360,
                -graph_title=> "",
                -x_title    => "",
                -y_title    => "",
                -startpos    => 0,
                -endpos        => -1,  # use $motifSize once it's known
                -greyscale    => 0,
                -regular_font => $REGFONT,
                -bold_font => $BOLDFONT,
                -title_font_size => 15,
                -axis_font_size  => 15,
                -height_logodds  => 0,
                @_);

    #read in passed arguments
    my ($xsize, $ysize, $xAxis, $yAxis, $title, $startpos, $endpos, $greyscale, $regular_font, $bold_font, $axis_font_size, $title_font_size, $height_logodds)
    = @args{qw(-xsize -ysize -x_title -y_title -graph_title -startpos -endpos -greyscale -regular_font -bold_font -axis_font_size -title_font_size -height_logodds)};

    my $data = $self->_raw_data($height_logodds);
    #initialize the most important variables first
    my $ICM = $data->{ICM};

    my $maxInfContent = $data->{maxIC};
    my $HPM = $data->{HPM}; # calculate Hitting Prob.
    my $width = $data->{width};

    my $motifSize = $ICM->getdim(0);
    my @alphabet = @{$self->alphabet};

    if ($endpos < 0 ) {
       $endpos = $motifSize;
    }


    #cut out the part of the profile that should be drawn - HPM first
    my $HPV = ones($motifSize);
    $HPV->slice('1:'.($motifSize-1).':2') .= $HPM->slice(':,1'); #get the original HPM-Values for Insert states
    #every position has 2 columns...
    $startpos *= 2;
    $endpos *= 2;
    #There should be at least 2 y-tics
    if($maxInfContent<2){
            $maxInfContent=2;
    }
    #check for start/endpos validity
    if($endpos>$motifSize){
        $endpos = $motifSize;
    }
    if($startpos>$endpos){
        die("Start position larger than end position!\n");
    }

    #cut out the part of the profile that should be drawn - Then the vectors
    $ICM = $ICM->slice("$startpos:".($endpos-1));
    $width = $width->slice("$startpos:".($endpos-1));
    $HPV = $HPV->slice("$startpos:".($endpos-1));

    #Imager-Objects
    my $black=NC(0,0,0);
    my $white=NC(255,255,255);
    #fill the bars for insert states with either a solid color or slashes
    my $hpFill = $greyscale?(Imager::Fill->new(hatch=>'dots16', fg=>$black, bg=>$white)):(Imager::Fill->new(solid => NC('#FF6666')));
    my $contribFill=$greyscale?(Imager::Fill->new(hatch=>'dots4', fg=>$black, bg=>$white)):(Imager::Fill->new(solid => NC('#FFCCCC')));

    #draw letters either with different colors or with 4 shades of grey
    my %letterColors;
    if($greyscale){
        if(@alphabet>4){
            %letterColors = ('A' => NC('#CCCCCC'),
                             'C' => NC('#CCCCCC'),
                             'D' => NC('#666666'),
                             'E' => NC('#666666'),
                             'F' => NC('#CCCCCC'),
                             'G' => NC('#333333'),
                             'H' => NC('#999999'),
                             'I' => NC('#CCCCCC'),
                             'K' => NC('#999999'),
                             'L' => NC('#CCCCCC'),
                             'M' => NC('#CCCCCC'),
                             'N' => NC('#333333'),
                             'P' => NC('#CCCCCC'),
                             'Q' => NC('#333333'),
                             'R' => NC('#999999'),
                             'S' => NC('#333333'),
                             'T' => NC('#333333'),
                             'V' => NC('#CCCCCC'),
                             'W' => NC('#CCCCCC'),
                             'Y' => NC('#333333')
                             );
        }
        else{
            %letterColors = ('A' => NC('#333333'),
                             'C' => NC('#666666'),
                             'G' => NC('#999999'),
                             'T' => NC('#CCCCCC'),
                             'U' => NC('#CCCCCC')
                             );
        }
    }
    else{
        if(@alphabet>4){
            %letterColors = ('A' => NC('#FF9966'),
                             'C' => NC('#009999'),
                             'D' => NC('#FF0000'),
                             'E' => NC('#CC0033'),
                             'F' => NC('#00FF00'),
                             'G' => NC('#FFFF00'),
                             'H' => NC('#660033'),
                             'I' => NC('#CC9933'),
                             'K' => NC('#663300'),
                             'L' => NC('#FF9933'),
                             'M' => NC('#CC99CC'),
                             'N' => NC('#336666'),
                             'P' => NC('#0099FF'),
                             'Q' => NC('#6666CC'),
                             'R' => NC('#990000'),
                             'S' => NC('#0000FF'),
                             'T' => NC('#00FFFF'),
                             'V' => NC('#FFCC33'),
                             'W' => NC('#66CC66'),
                             'Y' => NC('#006600')
                             );
        }
        else{
            %letterColors = ('A' => NC('#00FF00'),
                             'C' => NC('#0000FF'),
                             'G' => NC('#FFCC00'),
                             'T' => NC('#FF0000'),
                             'U' => NC('#FF0000')
                             );
        }
    }


    if ($args{'-letter_colours'} && ref($args{'-letter_colours'}) eq "HASH")
    {
        %letterColors = ();
        foreach my $key (keys %{$args{'-letter_colours'}})
        {
            $letterColors{$key} = NC($args{'-letter_colours'}->{$key})
        }
    }

    my $normFont=NF(file=>"$regular_font", type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";
    my $titleFont=NF(file=>"$bold_font", type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";
    my $vertFont=NF(file=>"$regular_font", type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";

    #make label vertical (transformation by rotation-matrix!)
    $vertFont->transform(matrix=>[ 0, -1, 0,
                                  1, 0, 0]);

    #make new image, fill white
    my $i=Imager->new(xsize=>$xsize, ysize=>$ysize, channels=>($greyscale?1:4)); # destination image
    $i->box(filled=>1,color=>$white);                     # fill with background color
                                                        # draw shadow just to right and below
                                                        # where log will be drawn

    #draw axis-descriptions and title
    $i->string(font=>$vertFont,
               string=>$yAxis,
               x=>18,
               y=>0.5*$ysize + 0.5*($vertFont->bounding_box(string=>$yAxis, size => $axis_font_size))[2],
               size=>$axis_font_size,
               color=>$black,
               aa=>1);

    $i->string(font=>$normFont,
               string=>$xAxis,
               x    =>  0.5*$xsize - 0.5*($normFont->bounding_box(string=>$xAxis, size => $axis_font_size))[2],
               y    =>  $ysize-7,
               size =>  $axis_font_size,
               color    =>  $black,
               aa=>1);

    $i->string(font=>$titleFont,
               string=>$title,
               x    =>  $xsize - ($titleFont->bounding_box(string=>$title, size => $title_font_size))[2]-10,
               y    =>  21,
               size =>   $title_font_size,
               color=>$black,
               aa=>1);

    #draw axes
    my $leftMargin = $yAxis?(($vertFont->bounding_box(string=>$yAxis))[3]-($vertFont->bounding_box(string=>$yAxis))[1]+20):20;
    my $rightMargin = $xsize - 24;
    my $lowerMargin = $xAxis?($ysize-($normFont->bounding_box(string=>$xAxis))[3]-($normFont->bounding_box(string=>$xAxis))[1]-25):($ysize-25);
    my $upperMargin = $title?(($titleFont->bounding_box(string=>$title))[3]-($titleFont->bounding_box(string=>$title))[1]+24):24;

    $i->line(color=>$black,
             x1=> $leftMargin,
             x2=> $rightMargin,
             y1=> $lowerMargin,
             y2=> $lowerMargin);

    $i->line(color=>$black,
             x1=> $leftMargin,
             x2=> $leftMargin,
             y1=> $upperMargin,
             y2=> $lowerMargin);

    #draw vertical tics and numbers
    my $yStep = ($lowerMargin-$upperMargin) / PDL::Math::ceil($maxInfContent);

    for(my $k=0; $k<=PDL::Math::ceil($maxInfContent); ++$k){
        $i->line(color=>$black,
                x1=>$leftMargin-3,
                x2=>$leftMargin+3,
                y1=>$lowerMargin - $k*$yStep,
                y2=>$lowerMargin - $k*$yStep);

        $i->string(font=>$normFont,
                   string=>"$k ",
                   x=>$leftMargin - 11,
                   y=>$lowerMargin - $k*$yStep +4,
                   size=>10,
                   color=>$black,
                   aa=>1);
    }

    #create Logo

    #calculate scaling-factor for given image-size
    my $scaleFactor = ($rightMargin - $leftMargin)/$width->sumover();
    #scale width to target-width
    $width *= $scaleFactor;
    #scale HP to target-width
    $HPV *= $scaleFactor;

    #print a Scale-Bar to show the relative sizes
    $i->box(color=>$black,
                xmin=>$leftMargin,
                xmax=>$leftMargin + $scaleFactor,
                ymin=>$upperMargin - 14,
                ymax=>$upperMargin - 15);
    $i->line(color=>$black,
                x1=>$leftMargin,
                x2=>$leftMargin,
                y1=>$upperMargin - 11,
                y2=>$upperMargin - 17);
    $i->line(color=>$black,
                x1=>$leftMargin + $scaleFactor,
                x2=>$leftMargin + $scaleFactor,
                y1=>$upperMargin - 11,
                y2=>$upperMargin - 17);

    my $isMatch = 1;

    for(my $k=0; $k<($endpos-$startpos); ++$k){
        my %infContent = ();
        #get all the values for one column from the matrix
        my @infValues = list $ICM->slice($k);

        my $l=0;

        #fill infContent-Hash with Values
        foreach my $key (@alphabet){
            if( $infValues[$l] !~ /^BAD$/ && ($infValues[$l] * $yStep) > 1){ # if the character would print lower than 1 pixel, discard it
                $infContent{$key} = sprintf('%.16f', $infValues[$l]);
            }
            ++$l; # this is the index of the information-content of the next character in the infValues vector
        }

        #define drawing order by descending infContent
        my @drawOrder = sort {$infContent{$a}<=>$infContent{$b}} keys(%infContent); # Sort the letters to be printed by decreasing inf-content

        #reset the height to 0
        my $sumHeight = $lowerMargin;

        unless($isMatch){
            #draw tick and number

            $i->line(color=>$black,
                    x1=>$leftMargin + $width->slice("0:$k")->sumover(),
                    x2=>$leftMargin + $width->slice("0:$k")->sumover(),
                    y1=>$lowerMargin - 3,
                    y2=>$lowerMargin + 3);
            $i->string(font=>$normFont,
                   string=>(($startpos + $k+1)/2).' ',
                   x=>$leftMargin + $width->slice("0:$k")->sumover() - 0.5*$width->slice(($k-1).':'.$k)->sumover() - 0.5*(($normFont->bounding_box(string=>($startpos + $k+1)/2))[2]-($normFont->bounding_box(string=>($startpos + $k+1)/2))[0]),
                   y=>$lowerMargin + 13,
                   size=>10,
                   color=>$black,
                   aa=>1);
            #draw bg for contribution
            $i->box(fill=>$contribFill, #use the fill selected earlier
                    xmin=>$leftMargin + $width->slice("0:".($k-1))->sumover(),
                    ymin=>$upperMargin,
                    xmax=>$leftMargin + $width->slice("0:$k")->sumover(),
                    ymax=>$lowerMargin-1);
            #draw bg for hitting probability
            $i->box(fill=>$hpFill, #use the fill selected earlier
                    xmin=>$leftMargin + $width->slice("0:".($k-1))->sumover(),
                    ymin=>$upperMargin,
                    xmax=>$leftMargin + $width->slice("0:".($k-1))->sumover() + $HPV->at($k),
                    ymax=>$lowerMargin-1);
            $isMatch = 1;
        }
        else{
            foreach my $key (@drawOrder){

                #get the size of the character to draw
                my @bbox = $normFont->bounding_box(string=>$key, size=>$infContent{$key}*$yStep, sizew=>$width->at($k,0));
                #calculate height and width of the character
                my $letterHeight = $bbox[5]-$bbox[4];
                my $letterWidth = $bbox[6]-$bbox[0];
                #there is a difference between target height and real drawn height, because many characters do not use all the space provided.
                #To solve that issue, we calculate a factor to multiply the target height with, so it compensates the difference in sizes.
                my $scaleFactor = $infContent{$key}* $yStep / $letterHeight;
                my $scaleFactorW = $width->at($k,0)/$letterWidth;
                if($infContent{$key}){
                    $i->string(font=>$normFont,
                           string=>$key,
                            # We want centered text; thus, we start from the middle of the column and go left by half of the positiv width of the letter.
                            # As letters tend to break out of their boundaries, we downscale them a little and only use 0.47 instead of 0.5
                           x=>$leftMargin + ($k?($width->slice("0:".($k-1))->sumover()):0) + (0.5*$width->at($k,0)-0.47*$bbox[6]*$scaleFactorW),
                           y=>$sumHeight + $bbox[4],
                           size=>$infContent{$key} * $yStep * $scaleFactor,
                           sizew=>$width->at($k,0)*$scaleFactorW*0.93,
                           color=>$letterColors{$key},
                           aa=>1);
                    $sumHeight -= $infContent{$key}*$yStep - ($bbox[4]*$scaleFactor);
                }
            }
            $isMatch = 0;
        }
    }
    $i->write(file=>$args{'-file'});
    return 1;
}


=head2 print_logo_dimensions

 Usage      : $pHMM->print_logo_dimensions(%args)
 Function   : print the dimensions used to build a logo from the represented HMM
 Returns    : true, if dimensions written properly, else false
 Args       :

    -file,          # name of the output file; should end in .png
    -xsize,         # Image width; default is 600
    -ysize          # Image heigth; default is 360
    -height_logodds  # by default, information content height is divided according to emission probabilities;
                     # if non-zero, heights will be based on log-odds scores (only positive scoring residues)

=cut

sub print_logo_dimensions
{
    my $self = shift;


    #default arguments
    my %args = (-xsize          => 600,
                -ysize          => 360,
                -height_logodds => 0,
                @_);

    #read in passed arguments
    my ($xsize, $ysize, $height_logodds)  = @args{qw(-xsize -ysize -height_logodds)};


    my $data = $self->_raw_data($height_logodds);
    #initialize the most important variables first
    my $ICM = $data->{ICM};

    my $maxInfContent = $data->{maxIC};
    my $HPM = $data->{HPM}; # calculate Hitting Prob.
    my $width = $data->{width};

    my $motifSize = $ICM->getdim(0);
    my @alphabet = @{$self->alphabet};


    #cut out the part of the profile that should be drawn - HPM first
    my $HPV = ones($motifSize);
    $HPV->slice('1:'.($motifSize-1).':2') .= $HPM->slice(':,1'); #get the original HPM-Values for Insert states

    #calculate scaling-factor for given image-size
    my $xScale = $xsize/$width->sumover();
    #scale width to target-width
    $width *= $xScale;
    #scale HP to target-width
    $HPV *= $xScale;


    my $yScale = $ysize / PDL::Math::ceil($maxInfContent);

    for(my $k=0; $k<$motifSize; $k += 2){

        my %infContent = ();
        #get all the values for one column from the matrix
        my @infValues = list $ICM->slice($k);
        my $l=0;

        #fill infContent-Hash with Values
        foreach my $key (@alphabet){
            $infContent{$key} = sprintf('%.6f', $infValues[$l]);
            ++$l; # this is the index of the information-content of the next character in the infValues vector
        }

        #define drawing order by descending infContent
        my @drawOrder = sort {$infContent{$a}<=>$infContent{$b}} keys(%infContent); # Sort the letters to be printed by decreasing inf-content

        printf "%d: insert_widths:%.3f,%.3f ; res_width:%.3f ; res_heights:",
                 ($k+2) / 2 ,
                 ( $width->slice("0:".($k+1))->sumover()) - ( $width->slice("0:".($k))->sumover() ) ,
                 ($width->slice("0:".($k))->sumover() + $HPV->at($k+1)) - ( $width->slice("0:".($k))->sumover()),
                 $width->at($k,0);

        foreach my $key (@drawOrder){
            if($infContent{$key} > 0 ){

                printf "$key=%.3f,", $infContent{$key}  * $yScale;

            }
        }

        print "\n";

    }


    return 1;

}
=head2 toHMMer2()

 Usage      : my $hmmer = $pHMM->toHMMer()
 Returns    : A string representing this HMM in HMMer file format

=cut

sub toHMMer2
{
    my $self = shift;

    my $file = 'HMMER'.$self->version."\n";
    $file .= "NAME  ".$self->name."\n";
    $file .= "ACC   ".$self->accession."\n" if $self->accession;
    $file .= "DESC  ".$self->description."\n" if $self->description;
    $file .= "LENG  ".$self->length."\n";
    $file .= "ALPH  ".(@{$self->alphabet}>4?'Amino':'Nucleic')."\n";
    $file .= "RF    no\n";
    $file .= "CS    no\n";
    $file .= "MAP   no\n";
    $file .= "NSEQ  ".$self->seqNumber."\n" if $self->seqNumber;
    $file .= "XT    ".HMM::Utilities::Score2Ascii($self->specialTransitionScores)." \n";
    $file .= "NULT ".HMM::Utilities::Score2Ascii($self->nullTransitionScores)."\n";
    $file .= "NULE ".HMM::Utilities::Score2Ascii($self->nullEmissionScores)." \n";
    $file .= "EVD   ".join('   ', list ($self->evidence) )."\n" if (defined($self->evidence) && dims($self->evidence));
    $file .= "HMM        ".join('      ', @{$self->alphabet})."    \n";
    $file .= "         m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e\n";
    $file .= '      '.HMM::Utilities::Score2Ascii($self->startTransitionScores)."\n";
    for (my $k=0; $k<$self->length; ++$k)
    {
        $file .= sprintf('%6d', $k+1);
        $file .= HMM::Utilities::Score2Ascii($self->emissionScores->slice("$k,:,0")).sprintf("%6d", $self->state2column($k+1)||($k+1))."\n"; # Match Emissions
        $file .= '     -'.HMM::Utilities::Score2Ascii($self->emissionScores->slice("$k,:,1"))." \n"; # Insert Emissions
        $file .= '     -'.HMM::Utilities::Score2Ascii($self->transitionScores->slice("$k,:"))." \n"; # Transitions
    }
    $file .= "//\n";
}

#private Function; Parse a file in HMMer-Format
sub _parseFile {
    my ($self, $file) = @_;

    # Boolean variables for proper parsing of *.hmm
    my $in_model = 0;
    my $got_trans = 0;
    my $in_state = 0;
    my $statenum = 0;    #Array Index
    foreach(split(/\n/, $file)) {
        if (/^NAME\s*(\S+)/) { $self->{'modelName'} = $1; }
        elsif (/^ACC\s*(\S+)/) { $self->{'accession'} = $1; }
        elsif (/^DESC\s*(.*)/) { $self->{'description'} = $1; }
        elsif (/^LENG\s*(\S+)/) { $self->{'length'} = $1; }
        elsif (/^ALPH\s*(\S+)/) { ($1 eq 'Amino') ?( $self->{'alphabet'}=20 ):( $self->{'alphabet'}=4 ) }
        elsif (/^NSEQ\s*(\S+).*/) { $self->{'seqNumber'} = $1; }
        elsif (/^XT\s*(.+)/) {
            $self->{'specialTransitions'} = HMM::Utilities::Ascii2Prob($1, 1);
            $self->{'specialTransitionScores'} = HMM::Utilities::Ascii2Score($1);
        }
        elsif (/^NULT\s+(\S.*)/) {
            $self->{'nullTransitions'} = HMM::Utilities::Ascii2Prob($1, 1);
            $self->{'nullTransitionScores'} = HMM::Utilities::Ascii2Score($1);
        }
        elsif (/^NULE\s+(\S.+)/) {
            $self->{'nullEmissions'} = HMM::Utilities::Ascii2Prob($1, 1/$self->alphabet);
            $self->{'nullEmissionScores'} = HMM::Utilities::Ascii2Score($1);
        }
        elsif (/^EVD\s*(\S+\s+\S+)/)
        {
            $self->{'evidence'} = HMM::Utilities::Ascii2Score($1);
        }
        elsif (/^\/\//) {
            last;
        }

        ############### start parsing the main-model ################

        elsif (/^HMM\s+([\w\s]*)/) {
            $in_model = 1; #we saw the first line of main model
            $self->{'alphabet'} = [split(/\s+/, $1)];
            #initialize data structures for HMM-representation
            $self->{'transitions'} = zeroes($self->{'length'}, 9); #A matrix holding the transition-probabilities
            $self->{'transitionScores'} = zeroes($self->{'length'}, 9); #A matrix holding transition-scores as given in the file
            $self->{'emissions'} = zeroes($self->{'length'}, scalar(@{$self->{'alphabet'}}), 2); #2 matrices holding the emission-prob. for M and I states
            $self->{'emissionScores'} = zeroes($self->{'length'}, scalar(@{$self->{'alphabet'}}), 2); #2 matrices holding the emission-scores for M and I states
        }

        elsif (!$got_trans && $in_model && (/\w->\w/)) {
            $got_trans = 1; #The human readable line which is the transition names
        }
        # The special B->D transition line
        elsif (!$in_state && $got_trans && (/\s*(.+)/)) {
            $self->{'startTransitions'} = HMM::Utilities::Ascii2Prob($1, 1.0);
            $self->{'startTransitionScores'} = HMM::Utilities::Ascii2Score($1);
            $in_state = 1;
        }
        elsif ($in_state && (/\s*(\S+)\s*(.+)/)) {
            my ($s, $l) = ($1, $2);
            my @split = split(/\s+/, $l);
            if ($s =~ /\d+/) { # Match Line
                $statenum = int($s)-1;
                my @res = $l =~ /([\-\d]+)\s*/g;
                # If there's a column that maps the state position to the alignment, save it
                if (@res > @{$self->{'alphabet'}})
                {
                    my $col = pop @res;
                    $self->{'state2column'} ||= [];
                    $self->{'state2column'}->[$statenum] = $col;
                }
                $self->{'emissions'}->slice("$statenum,:,0") .= HMM::Utilities::Ascii2Prob(join (' ', @res), $self->{'nullEmissions'});
                $self->{'emissionScores'}->slice("$statenum,:,0") .= HMM::Utilities::Ascii2Score(join (' ', @res));
            }
            elsif (scalar(@split) == 9) { # transition line
                $self->{'transitions'}->slice("$statenum") .= HMM::Utilities::Ascii2Prob($l, 1.0);
                $self->{'transitionScores'}->slice("$statenum") .= HMM::Utilities::Ascii2Score($l);
            }
            else { # Insert Line
                $self->{'emissions'}->slice("$statenum,:,1") .= HMM::Utilities::Ascii2Prob($l, $self->{'nullEmissions'});
                $self->{'emissionScores'}->slice("$statenum,:,1") .= HMM::Utilities::Ascii2Score($l);
            }
        }
    }
    #
    # Check that the length is correct
    #
    if ($self->length() != $self->emissions->getdim(0))
    {
        warn("Length given in annotation is not identical to number of states in model, correcting...");
        $self->{'length'} = $self->emissions->getdim(0);
    }
    return 1;
}

sub _parseFile3 {
    my ($self, $file) = @_;

    my @lines = split(/\n/, $file);
    shift @lines; # dump first line which contains the HMMER version description
    my $line;
    do {
        $line = shift @lines;
        if ($line =~ /^NAME\s*(\S+)/) { $self->{'modelName'} = $1; }
        elsif ($line =~ /^ACC\s*(\S+)/) { $self->{'accession'} = $1; }
        elsif ($line =~ /^DESC\s*(.*)/) { $self->{'description'} = $1; }
        elsif ($line =~ /^LENG\s*(\S+)/) { $self->{'length'} = $1; }
        elsif ($line =~ /^ALPH\s*(\S+)/) { ($1 eq 'Amino') ?( $self->{'alphabet'}=20 ):( $self->{'alphabet'}=4 ) }
        elsif ($line =~ /^NSEQ\s*(\S+).*/) { $self->{'seqNumber'} = $1; }
        #
        # There are a lot of new HMMER3 specific tags missing here
        #
    } while (defined $line and $line !~ /^HMM/);

    unless ($self->length and (defined $self->name or defined $self->accession))
    {
        warn("Parse error before $line");
        return;
    }

    if ($line =~ /^HMM\s+([\w\s]*)/)
    {
        $self->{'alphabet'} = [split(/\s+/, $1)];
        #cleanly initialize data structures for HMM-representation
        $self->{'transitions'} = zeroes($self->length, 7); #A matrix holding the transition-probabilities
        $self->{'transitionScores'} = zeroes($self->length, 7); #A matrix holding transition-scores as given in the file
        $self->{'emissions'} = zeroes($self->length, scalar(@{$self->{'alphabet'}}), 2); #2 matrices holding the emission-prob. for M and I states
        $self->{'emissionScores'} = zeroes($self->length, scalar(@{$self->{'alphabet'}}), 2); #2 matrices holding the emission-scores for M and I states
    }
    else
    {
        warn("HMM is badly formatted");
        return;
    }
    $line = shift @lines;
    unless ($line =~ /m->m/)
    {
        warn("HMM is badly formatted");
        return;
    }
    $line = shift @lines;
    # Optional COMPO line
    if ($line =~ /^\s*COMPO\s+(\S.*)/)
    {
        $self->{'averageComposition'} = HMM::Utilities3::Ascii2Prob($1)->transpose();
        $line = shift @lines;
    }
    # Parse the null model
    if ($line =~ /^\s*(\S.*\S)\s*$/)
    {
        $self->{'nullEmissions'} = HMM::Utilities3::Ascii2Prob($1)->transpose();
        $self->{'nullEmissionScores'} = HMM::Utilities3::Ascii2Log($1)->transpose();
    }
    else
    {
        warn("HMM is badly formatted");
        return;
    }
    # Parse the initial transitions to the first match state
    $line = shift @lines;
    if ($line =~ /^\s*(\S.*\S)\s*$/)
    {
        $self->{'startTransitions'} = HMM::Utilities3::Ascii2Prob($1)->transpose();
        $self->{'startTransitionScores'} = HMM::Utilities3::Ascii2Log($1)->transpose();
    }
    else
    {
        warn("HMM is badly formatted");
        return;
    }

    # Now go through the model section and parse each state, consisting of 3 lines ()
    my $substate = 0;
    my $stateno;
    $line = shift @lines;
    while (defined $line and $line !~ m|//|)
    {
        # Match line
        if ($substate == 0)
        {
            my ($s, $ascii, $map, $consensus, $rf, $cs);

            if ($self->version =~ /[abcd]$/) {
               ($s, $ascii, $map, $rf, $cs)             = $line =~ /^\s*(\d+)\s+(.*)\s+(\S+)\s+(\S)\s+(\S)\s*$/;
            } else {  # hmm version e or greater
               ($s, $ascii, $map, $consensus, $rf, $cs) = $line =~ /^\s*(\d+)\s+(.*)\s+(\S+)\s+(\S)\s+(\S)\s+(\S)\s*$/;
            }



            unless (defined $s and defined $map and defined $rf and defined $cs and defined $ascii)
            {
                warn("Parse error: $line");
                return;
            }
            $stateno = $s-1;
            $self->{state2column}[$stateno] = $map unless $map eq '-';
            $self->{consensus} .= $cs unless $cs eq '-';
            $self->{reference} .= $rf unless $rf eq '-';

            $self->{'emissions'}->slice("$stateno,:,0") .= HMM::Utilities3::Ascii2Prob($ascii)->transpose();
            $self->{'emissionScores'}->slice("$stateno,:,0") .= HMM::Utilities3::Ascii2Log($ascii)->transpose();
            $substate++;
        }
        # Insert line
        elsif ($substate == 1)
        {
            my ($ascii) = $line =~ /^\s*(\S.*\S)\s*$/;
            unless (defined $ascii)
            {
                warn("Parse error: $line");
                return;
            }
            $self->{'emissions'}->slice("$stateno,:,1") .= HMM::Utilities3::Ascii2Prob($ascii)->transpose();
            $self->{'emissionScores'}->slice("$stateno,:,1") .= HMM::Utilities3::Ascii2Log($ascii)->transpose();
            $substate++;
        }
        # Transition line
        else
        {
            my ($ascii) = $line =~ /^\s*(\S.*\S)\s*$/;
            $self->{'transitions'}->slice("$stateno,:") .= HMM::Utilities3::Ascii2Prob($ascii)->transpose();
            $self->{'transitionScores'}->slice("$stateno,:") .= HMM::Utilities3::Ascii2Log($ascii)->transpose();
            $substate = 0;
            $stateno = undef;
        }
        $line = shift @lines;
    }

    #
    # Check that the length is correct
    #
    if ($self->length() != $self->emissions->getdim(0))
    {
        warn("Length given in annotation is not identical to number of states in model, correcting...");
        $self->{'length'} = $self->emissions->getdim(0);
    }
    return 1;
}

return 1;

=head1 BUGS

This is an alpha release; if you witness malfunction please contact the author!

=head1 AUTHOR

Benjamin Schuster-Boeckler <boeckler(at)molgen.mpg.de>

=head1 SEE ALSO

L<HMM>, L<HMM::Utilities>, L<TFBS::Matrix>, L<Imager>, L<PDL>

=cut
