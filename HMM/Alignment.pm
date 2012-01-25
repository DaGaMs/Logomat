# This is Clan::Alignment written by Benjamin Schuster-Boeckler
# Released under the same license as Perl itself


=head1 NAME

HMM - representation of hidden markov moddels for computational biology

=head1 SYNOPSIS

    use HMM::Alignment;
    
    my $alignment = HMM::Alignment->new();
    
    $alignment->align(-HMM1id   => 'FGF',
                      -HMM2id   => 'Agglutinin');
    
    $alignment->draw_logo(-file  => 'fgf-aggl_al.png',
                             -xsize => 20*100,
                             -ysize => 400,
HMM2                             -y_title1 => 'FGF',
                             -y_title2 => 'Aggl',
                             -graph_title => 'FGF vs. Aggl') || die "Failed drawing logo: $!\n";


=head1 DESCRIPTION

This class serves as a wrapper around PRC, the Profile Comparer, written by Martin Madera (http://supfam.mrc-lmb.cam.ac.uk/PRC/). 
It can read a file in PRC format or run PRC directly, if it is installed on the system. It needs a pair of HMM objects to do so.
The user can either provide these objects or let the module create them from HMM files or PFAM ids. There is a method that creates
a graphical representation of the alignment. It is an extension of the HMM Logo concept and follows the same concepts.

=head1 METHODS

=cut

package HMM::Alignment;

use strict;
use warnings;
use vars '@ISA'; #inheritance
use PDL::LiteF;
use Imager ':handy';
use Imager::Fill;
use HMM::Profile;
use HMM::Utilities;
use IO::File;
use POSIX qw(:sys_wait_h tmpnam); # We need this for temporary file storage!
use Env qw(PRC_BIN PATH);

our $REGFONT = '/usr/local/share/fonts/ttfonts/arial.ttf';
our $BOLDFONT = '/usr/local/share/fonts/ttfonts/arialbd.ttf';

=head2 new(%args)

 Usage      : my $alignment = HMM::Alignment->new(%args)
 Function : constructor for the HMM::Alignment object
 Returns  : a new HMM::Alignment object
 Args      :
    -HMM1file    # filename or -handle, HMMer plan7 format
    -HMM2file    # filename or -handle, HMMer plan7 format
    -HMM1id    # PFAM id; The module will try to fetch the HMM from PFAM through LWP
    -HMM2id    # PFAM id
    -HMM1       # HMM object
    -HMM2       # HMM object
    -prcfile    # Output of PRC in "prc" format

=cut

# Create a new HMM-Aligment Object. This is merely an abstraction of the output from PRC, at the moment.
# Additionally, we need 2 HMM objects.
sub new
{
    my ($class, %args) = @_;
    my $self = {};
    bless($self, ref($class) || $class);
    if(ref($args{'-prcfile'}) eq 'Fh' or ref($args{'-prcfile'}) eq 'GLOB' or ref($args{'-prcfile'}) eq 'IO::File')
    {
        my $file;
        my $fh = $args{'-prcfile'};
        while(<$fh>)
        {
            $file .= $_;
        }
        $self->_parseFile($file);
    }
    elsif ($args{'-prcfile'})
    {
        my $file;
        open(HMM, $args{'-prcfile'}) || (warn("Couldn't open prc file: $!\n") && return 0);
        while(<HMM>){
            $file .= $_;
        }
        close(HMM);
        $self->_parseFile($file);
    }
    
    # Try to create 2 HMM objects. First, try the name from the alignment files. Then, check the passed arguments. If no HMM is there, die!
    $self->{'HMM1'} = HMM::Profile->new( '-hmmerfile' => $args{'-HMM1file'} ) if $args{'-HMM1file'};
    $self->{'HMM1'} = HMM::Profile->new( '-pfamid' => $args{'-HMM1id'} ) if $args{'-HMM1id'};
    $self->{'HMM1'} = $args{'-HMM1'} if $args{'-HMM1'};
    $self->{'HMM1'} = HMM::Profile->new( '-pfamid' => $self->{'hmm1'} ) unless $self->{'HMM1'};
    (warn("Missing HMM in arguments!\n") && return 0) unless $self->{'HMM1'};
    
    $self->{'HMM2'} = HMM::Profile->new( '-hmmerfile' => $args{'-HMM2file'} ) if $args{'-HMM2file'};
    $self->{'HMM2'} = HMM::Profile->new( '-pfamid' => $args{'-HMM2id'} ) if $args{'-HMM2id'};
    $self->{'HMM2'} = $args{'-HMM2'} if $args{'-HMM2'};
    $self->{'HMM2'} = HMM::Profile->new( '-pfamid' => $self->{'hmm2'} ) unless $self->{'HMM2'};
    (warn("Missing HMM in arguments!\n") && return 0) unless $self->{'HMM2'};
    
    if( $self->{'HMM1'} and $self->{'HMM2'} and $self->{'start1'} and $self->{'end1'} and $self->{'start2'} and $self->{'end2'} )
    {
        my $startPos1 = $self->{'start1'} - 1; #First position has index 0, but PRC gives it index 1
        my $endPos1 = $self->{'end1'} - 1;
        my $startPos2 = $self->{'start2'} - 1;
        my $endPos2 = $self->{'end2'} - 1;
        $self->{'HMM1'}->{'emissions'} = $self->{'HMM1'}->emissions->slice("$startPos1:$endPos1,:,:");
        $self->{'HMM2'}->{'emissions'} = $self->{'HMM2'}->emissions->slice("$startPos2:$endPos2,:,:");
        $self->{'HMM1'}->{'transitions'} = $self->{'HMM1'}->transitions()->slice("$startPos1:$endPos1,:");
        $self->{'HMM2'}->{'transitions'} = $self->{'HMM2'}->transitions()->slice("$startPos2:$endPos2,:");
    }

    return $self;
}

=head2 align(%args)

 Usage      : my $alignment = HMM::Alignment->align(%args)
 Function : This method either aligns the HMMs given at construction of the object (if called /out arguments) or
            creates a new object and aligns the HMMs passed as arguments
 Returns  : an HMM::Alignment object
 Args      :
    -HMM1file    # filename or -handle, HMMer plan7 format
    -HMM2file    # filename or -handle, HMMer plan7 format
    -HMM1id    # PFAM id; The module will try to fetch the HMM from PFAM through LWP
    -HMM2id    # PFAM id
    -HMM1       # HMM object
    -HMM2       # HMM object
    -algorithm  # Algorithm to use; Either 'forward' (def) or 'viterbi'
    -mode       # 'local-local' (def), 'global-global' or 'local-global' (and vice-versa)
    -scoreFunc  # 'dot1' (def) or dot2

=cut

# This method takes 2 HMMs (or the ones that the object was initialized with) and aligns them using prc
sub align {
    my ($self, %args) = @_;
    
        
    # Called as constructor for new object whenever this is a) called as a constructor or b) with explicit arguments
    if ( !ref($self) and values(%args) or values(%args) )
    {
        # Call real constructor
        unless($self = new($self, %args))
        {
            warn("Creation of alignment object failed!"); 
            return 0;
        }
    } 
    else
    {
        die("No HMM1 given") unless $self->{'HMM1'};
        die("No HMM2 given") unless $self->{'HMM2'};
    }

    # If there were no args given, then we should align the content of this object. 
    # If we're given a Pfam ID or a Filehandle, we have to write these to a file to align
    my ($HMM1, $HMM2);
    # try new temporary filenames until we get one that didn't already exist
    my ($fh1, $fh2);
    do
    {
        $HMM1 = tmpnam().".$$.tmp.hmm";
    } until $fh1 = IO::File->new($HMM1, O_RDWR|O_CREAT|O_EXCL);
    do
    {
        $HMM2 = tmpnam().".$$.tmp.hmm";
    } until $fh2 = IO::File->new($HMM2, O_RDWR|O_CREAT|O_EXCL);
    
    # install atexit-style handler so that when we exit or die,
    # we automatically delete this temporary file
    END
    { 
        unlink($HMM1) or die "Couldn't unlink $HMM1: $!\n" unless (!$HMM1 || !-e $HMM1 || $HMM1 !~ /$$.tmp.hmm/);
        unlink($HMM2) or die "Couldn't unlink $HMM2: $!\n" unless (!$HMM2 || !-e $HMM2 || $HMM2 !~ /$$.tmp.hmm/);
    }
    
    print $fh1 $self->{'HMM1'}->toHMMer();
    $fh1->close;
    print $fh2 $self->{'HMM2'}->toHMMer();
    $fh2->close;
    
    # Now run prc, if it is in the path        
    # Detaint the imported variables first
    if ($PRC_BIN)
    {
        ($PRC_BIN) = $PRC_BIN =~ /([a-z0-9\/\.\-_]+)/i;
    }
    else
    {
        $PRC_BIN = 'prc';
    }
    
    ($PATH) = $PATH =~ /([a-z0-9\/\.\-_:]+)/i if $PATH;
    
    if ($PATH && ! -x $PRC_BIN)
    {
        my @path = split (':', $PATH);
        do
        {
            $PATH = shift(@path);
            ($PATH) = $PATH =~ /(.*)\/?/ if $PATH; # remove trailing /
        } while ($PATH && ! -x "$PATH/$PRC_BIN" );
        $PRC_BIN = "$PATH/$PRC_BIN" if -x "$PATH/$PRC_BIN"; 
    }
    
    %args = (
        '-algorithm'    => 'viterbi',
        '-mode'         => 'local-local',
        '-scoreFunc'    => 'dot2',
        %args
    );
    
    my ($algorithm, $mode, $scoreFunc)
    = @args{qw(-algorithm -mode -scoreFunc)};
    
    # We have seen that on some machines, you need to do 2>&1 for some strange reason... Maybe a bug in Apache?
    open (IN, "$PRC_BIN -align prc -hits 1 -algo $algorithm -mode $mode -MMfn $scoreFunc  $HMM1 $HMM2 |") or die "Error: $!\n";
    
    my $string;
    while(<IN>)
    {
        $string .= $_;
    }
    close(IN);
    
    (warn("PRC returned empty result!") && return 0) unless ($string);

    $self->_parseFile($string);

    (warn("Missing values in PRC output or no overlap!") && return 0) unless ($self->{'string1'} and $self->{'string2'} and $self->{'start1'} and $self->{'start2'} and $self->{'end1'} and $self->{'end2'});

    my $startPos1 = $self->{'start1'} - 1; #First position has index 0, but PRC gives it index 1
    my $endPos1 = $self->{'end1'} - 1;
    my $startPos2 = $self->{'start2'} - 1;
    my $endPos2 = $self->{'end2'} - 1;
    $self->{'HMM1'}->{'emissions'} = $self->{'HMM1'}->emissions->slice("$startPos1:$endPos1,:,:");
    $self->{'HMM2'}->{'emissions'} = $self->{'HMM2'}->emissions->slice("$startPos2:$endPos2,:,:");
    $self->{'HMM1'}->{'transitions'} = $self->{'HMM1'}->transitions()->slice("$startPos1:$endPos1,:");
    $self->{'HMM2'}->{'transitions'} = $self->{'HMM2'}->transitions()->slice("$startPos2:$endPos2,:");

    return $self;
}

=head2 draw_logo()

 Usage      : my $alignment->draw_logo(%args)
 Function : Creates an image file visualizing the alignment in Alignment-Logo format
 Args     :
    -file           # Required: A filename to write the image to
    -xsize          # Horizontal width of the image
    -ysize          # Vertical height of the image
    -x_title        # X-Axis description
    -y_title1       # Y-Axis description of upper logo
    -y_title2       # Y-Axis description of lower logo
    -graph_title    # Title of the image (Will appear in upper right corner)
    -greyscale      # Not yet operational

=cut

# This is mainly a copy of drawLogo from HMM::Profile. The difference is that 2 logos are drawn
# simultaneously. 
sub draw_logo {
    my $self = shift;
    my $HMM1 = $self->{'HMM1'}; # Make a copy of the piddles and slice away the unnecessary stuff...
    my $HMM2 = $self->{'HMM2'};
    my $startPos1 = $self->{'start1'} - 1; #First position has index 0, but PRC gives it index 1
    my $endPos1 = $self->{'end1'} - 1;
    my $startPos2 = $self->{'start2'} - 1;
    my $endPos2 = $self->{'end2'} - 1;
    
    my $alignmentLength = _max(length($self->{'string1'}), length($self->{'string2'}));
    
    #default arguments
    my %args = (-ysize      => 600,
                -xsize        => 25*$alignmentLength,    
                -graph_title=> "",
                -x_title    => "",
                -y_title1    => "",
                -y_title2    => "",
                -greyscale    => 0,
                @_);

    #read in passed arguments
    my ($xSize, $ySize, $xAxis, $yAxis1, $yAxis2, $title, $greyscale)
    = @args{qw(-xsize -ysize -x_title -y_title1 -y_title2 -graph_title -greyscale)};

    #Imager-Objects
    my $black=NC(0,0,0);
    my $white=NC(255,255,255);

    my $normFont=NF(file=>$REGFONT, type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";
    my $titleFont=NF(file=>$BOLDFONT, type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";
    my $vertFont=NF(file=>$REGFONT, type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";

    my $normFontSize = 15; # Should this be dynamical?
    my $smallFontSize = 12;
    my $tinyFontSize = 10;
    
    #make label vertical (transformation by rotation-matrix!)
    $vertFont->transform(matrix=>[ 0, -1, 0, 
                                  1, 0, 0]);

    #make new image, fill white
    my $i=Imager->new(xsize=>$xSize, ysize=>$ySize, channels=>($greyscale?1:4)); # destination image
    $i->box(filled=>1,color=>$white);                     # fill with background color


    # Set the margins of the used image area
    my $leftMargin = ($yAxis1 || $yAxis2)?($vertFont->bounding_box(string=>$yAxis1, size=>$smallFontSize)->global_ascent - $vertFont->bounding_box(string=>$yAxis1, size=>$smallFontSize)->global_descent +20):20;
    my $rightMargin = $xSize - 20;
    my $topMargin = $title?($titleFont->bounding_box(string=>$title, size=>$normFontSize)->global_ascent - $titleFont->bounding_box(string=>$title, size=>$normFontSize)->global_descent +21):20;
    my $bottomMargin = $xAxis?($ySize- $normFont->bounding_box(string=>$xAxis, size=>$normFontSize)->global_ascent - $normFont->bounding_box(string=>$xAxis, size=>$normFontSize)->global_descent -25):($ySize-22);
    my $logoHeight = ($bottomMargin - $topMargin - 0.1*$ySize )/2; # The height of a single logo is a half of the space dedicated to one logo minus the space needed for the gap between the 2
    my $middleMargin = $bottomMargin - $logoHeight;
    
    #print the title
    $i->string(font=>$titleFont,
               string=>$title,
               x=>$xSize - $titleFont->bounding_box(string=>$title)->end_offset - 18, 
               y=>18,
               size=>$normFontSize,
               color=>$black,
               aa=>1);
               
    #Print the Y-Axis descriptions
    $i->string(font=>$vertFont,
               string=>$yAxis2,
               x=>18, # Should this be dynamical?
               y=>$middleMargin + 0.5 * $logoHeight + 0.5* $vertFont->bounding_box(string=>$yAxis2, size=>$smallFontSize)->end_offset,
               size=>$smallFontSize,
               color=>$black,
               aa=>1);
    $i->string(font=>$vertFont,
               string=>$yAxis1,
               x=>18, # Should this be dynamical?
               y=>$middleMargin - 0.5 * $logoHeight + 0.5* $vertFont->bounding_box(string=>$yAxis1, size=>$smallFontSize)->end_offset,
               size=>$smallFontSize,
               color=>$black,
               aa=>1);

    # Print the X-Axis description
    $i->string(font=>$normFont,
               string=>$xAxis,
               x=>0.5*$xSize - 0.5 * $normFont->bounding_box(string=>$xAxis, size=>$normFontSize)->end_offset, 
               y=>$ySize-3,
               size=>$normFontSize,
               color=>$black,
               aa=>1);
    
    ########## Start generating the logos themselves #############
    #initialize the most important variables first
    my $ICM1 = toICM($HMM1->emissions(), $HMM1->nullEmissions()); # Calculate Information Content Matrix
    my $ICM2 = toICM($HMM2->emissions(), $HMM2->nullEmissions()); # Calculate Information Content Matrix
    my $motifSize1 = $ICM1->getdim(0);
    my $motifSize2 = $ICM2->getdim(0);    

    # Calculate Hitting Probability
    my $HPM1 = toHPM($HMM1->startTransitions(), $HMM1->transitions()); 
    my $HPM2 = toHPM($HMM2->startTransitions(), $HMM2->transitions());

    # scale the Insert states relative to the estimated retention time, which is 1/(1-p(I->I)) = 1/p(I->M)
    my $width1 = zeroes($motifSize1);
    $width1->slice("0:".($motifSize1-2).":2") .= $HPM1->slice(':,0');
    $width1->slice("1:".($motifSize1-1).":2") .= ($HPM1->slice(':,1')/$HMM1->transitions()->slice(':,3'))->badmask(0); #scale Insert-width by estimated retention period
    my $width2 = zeroes($motifSize2);
    $width2->slice("0:".($motifSize2-2).":2") .= $HPM2->slice(':,0');
    $width2->slice("1:".($motifSize2-1).":2") .= ($HPM2->slice(':,1')/$HMM2->transitions()->slice(':,3'))->badmask(0); #scale Insert-width by estimated retention period
    
    #calculate scaling-factor for given image-size. Only consider the relevant states
    my $scaleFactor = ($rightMargin - $leftMargin)/ _max($width1->sumover(), $width2->sumover());     

    #print a Scale-Bar to show the relative sizes
    $i->box(color=>$black,
                xmin=>$leftMargin,
                xmax=>$leftMargin + $scaleFactor,
                ymin=>$topMargin - 14,
                ymax=>$topMargin - 15);
    $i->line(color=>$black,
                x1=>$leftMargin,
                x2=>$leftMargin,
                y1=>$topMargin - 11,
                y2=>$topMargin - 17);
    $i->line(color=>$black,
                x1=>$leftMargin + $scaleFactor,
                x2=>$leftMargin + $scaleFactor,
                y1=>$topMargin - 11,
                y2=>$topMargin - 17);
    
    # Now, create the 2 logos and paste them into the picture at the correct positions...
    $i->paste(left  => $leftMargin-15,
              top   => $topMargin,
              img   => $self->_createLogo($ICM1, $HPM1, $width1, $HMM1->alphabet, $startPos1*2, 20 + $width1->sumover() * $scaleFactor, 20 + $logoHeight, $self->{'string1'}, $greyscale, 0)
              );
              
    $i->paste(left  => $leftMargin-15,
              top   => $middleMargin,
              img   => $self->_createLogo($ICM2, $HPM2, $width2, $HMM2->alphabet, $startPos2*2, 20 + $width2->sumover() * $scaleFactor, 20 + $logoHeight, $self->{'string2'}, $greyscale, 1)
              );

    # Now we have to draw the connecting lines between aligned states
    my $stateString1 = $self->{'string1'};
    my $stateString2 = $self->{'string2'};
    $stateString1 =~ s/[D~]//g;
    $stateString2 =~ s/[D~]//g;
    my $lineCount1 = length($stateString1);
    my $lineCount2 = length($stateString2);
    die "Error in alignment: Unequal number of matched states!\n" unless $lineCount1 == $lineCount2; # Assertion
    $stateString1 = $self->{'string1'};
    $stateString2 = $self->{'string2'};
    $stateString1 =~ s/~//g; # Ignore gaps, they don't change anything about the connecting lines
    $stateString2 =~ s/~//g;
    
    #Create an array that tells us which states to pair up
    my $alignedStates = zeroes($lineCount1, 2, 2); # We need a third dimension for the x_from/x_to positions
    my $curState = 0;
    my $curPos = 0;

    #$curState += length($1) * 2 - 2 if ($stateString1 =~ s/^(D+)//);
    
    # Deal with pos 0
    if(unpack('A1', $stateString1) eq 'M') # We start with a match (normal case)
    {
        $alignedStates->slice("0,0,1") .= $width1->at(0);
    }
    elsif(unpack('A1', $stateString1) eq 'I') # We start with an insert (rare)
    {
        $alignedStates->slice("0,0,0") .= $width1->at(0);
        $alignedStates->slice("0,0,1") .= $width1->at(0) + $width1->at(1);
    }
    else # We start with a delete (insane)
    {
        $curPos = -1;
    }    

    for (my $k = 1; $k<length($stateString1); ++$k)
    {
        my ($lastChar, $alignedChar) = unpack('x'.($k-1).' A1 A1', $stateString1);
        
        if ($alignedChar eq 'M')
        {
            $curPos++; # Which position in the paired states array are we in?
            if ($lastChar eq 'I')
            {
                $curState++; # Which state in the logo are we dealing with?
            }
            else
            {
                $curState += 2;
            }
            eval{
            # Now that we have determined which states align, insert coordinates into matrix
            $alignedStates->slice("$curPos,0,0") .= $width1->slice("0:".($curState-1))->sumover();
            $alignedStates->slice("$curPos,0,1") .= $width1->slice("0:$curState")->sumover();
            };
            if($@)
            {
                warn("FATAL: State: $curState, Pos: $curPos, k: $k, Msg: ", $@);
                #return(0);
            }
        }
        elsif ($alignedChar eq 'I')
        {
            $curPos++;
            if ($lastChar ne 'I')
            {
                $curState++;
            }
            eval{
            # Now that we have determined which state align, isert coordinates into matrix
            $alignedStates->slice("$curPos,0,0") .= $width1->slice("0:".($curState-1))->sumover();
            $alignedStates->slice("$curPos,0,1") .= $width1->slice("0:$curState")->sumover();
            };
            if($@)
            {
                warn("FATAL: State: $curState, Pos: $curPos, k: $k, Msg: ", $@);
                #return(0);
            }
        }
        elsif ($alignedChar eq 'D')
        {
            $curState += 2;
        }
    }

    # Now do it all again for the second HMM
    $curState = 0;
    $curPos = 0;
    
    # Deal with pos 0
    if(unpack('A1', $stateString2) eq 'M')
    {
        $alignedStates->slice("0,1,1") .= $width2->at(0); # We start with a match (normal case)
    }
    elsif(unpack('A1', $stateString2) eq 'I')
    {
        $alignedStates->slice("0,1,0") .= $width2->at(0);
        $alignedStates->slice("0,1,1") .= $width2->at(0) + $width2->at(1); # We start with an insert (rare)
    }
    else # We start with a Gap/Delete (insane)
    {
        $curState++;
    }


    for (my $k = 1; $k<length($stateString2); ++$k)
    {
        my ($lastChar, $alignedChar) = unpack('x'.($k-1).' A1 A1', $stateString2);
        
        if ($alignedChar eq 'M')
        {
            $curPos++; # Which position in the paired states array are we in?
            if ($lastChar eq 'I')
            {
                $curState++; # Which state in the logo are we dealing with?
            }
            else
            {
                $curState += 2;
            }
            eval{
            # Now that we have determined which state align, isert coordinates into matrix
            $alignedStates->slice("$curPos,1,0") .= $width2->slice("0:".($curState-1))->sumover();
            $alignedStates->slice("$curPos,1,1") .= $width2->slice("0:$curState")->sumover();
            };
            if($@)
            {
                warn("FATAL: $@");
                #return(0);
            }
        }
        elsif ($alignedChar eq 'I')
        {
            $curPos++;
            if ($lastChar ne 'I')
            {
                $curState++;
            }
            eval{
            # Now that we have determined which state align, isert coordinates into matrix
            $alignedStates->slice("$curPos,1,0") .= $width2->slice("0:".($curState-1))->sumover();
            $alignedStates->slice("$curPos,1,1") .= $width2->slice("0:$curState")->sumover();
            };
            if($@)
            {
                warn("FATAL: $@");
                #return(0);
            }
        }
        elsif ($alignedChar eq 'D')
        {
            $curState += 2;
        }
    }
 
    $alignedStates *= $scaleFactor;
    $alignedStates += $leftMargin;

    my $whiteFill = Imager::Fill->new(solid=>$white, combine=>'normal');

    # Now, go through the matrix and draw the lines
    for (my $k = 0; $k < $alignedStates->getdim(0); ++$k)
    {
        $i->polygon(color   => $black,
                    'x' => [$alignedStates->at($k,0,0)+1, $alignedStates->at($k,0,1)-1, $alignedStates->at($k,1,1)-1, $alignedStates->at($k,1,0)+1],
                    'y' => [$topMargin + $logoHeight + 15, $topMargin + $logoHeight + 15, $middleMargin+5, $middleMargin+5],
                    aa  => 0
                    );
#         $i->line(color=>$black,
#                     x1  => $alignedStates->at($k,0,0),
#                     x2  => $alignedStates->at($k,1,0),
#                     y1  => $topMargin + $logoHeight + 15,
#                     y2  => $middleMargin+5,
#                     aa  => 0);
#         $i->line(color=>$black,
#                     x1  => $alignedStates->at($k,0,1),
#                     x2  => $alignedStates->at($k,1,1),
#                     y1  => $topMargin + $logoHeight + 15,
#                     y2  => $middleMargin+5,
#                     aa  => 0);
    }

    # Finally, draw the logo
    if ($args{'-file'}) {
        $i->write(type=>'png', file=>$args{'-file'});
    }
    elsif ($args{'-fh'}) {
        $i->write(type=>'png', fh=>$args{'-file'});
    }
    elsif ($args{'-data'}) {
        $i->write(type=>'png', data=>$args{'-data'});
    }
    else {
        warn __PACKAGE__." WARN No data output handle passed!\n";
        return;
    }
}

=head2 toPRC

  Usage     : $string = $alignment->toPRC()
  Function  : Returns a string representing this Alignment in PRC format

=cut

sub toPRC
{
    my $self = shift;
    my $out = "PRC ".$self->version." run through HMM::Alignment ".$HMM::VERSION.
              "\nCopyright (C) 2002-4 Martin Madera and MRC LMB, Cambridge, UK\nFreely distributed under the GNU General Public License\n";
    $out .= "Command     : ".$self->{'command'}."\n" if $self->{'command'};
    $out .= "Algorithm   : ".$self->algorithm."\n";
    $out .= "Align mode  : ".$self->mode."\n";
    $out .= "Alignments  : prc\n";
    $out .= "Simple Stop : ".$self->{'stop'}."\n" if $self->{'stop'};
    $out .= "Max hits    : 1\n";
    $out .= "Started     : ".$self->{'date'}."\n" if $self->{'date'};
    $out .= "\n#  hmm1  start1    end1 length1    hit_no    hmm2  start2    end2 length2  co-emis  simple  reverse\n";
    
    my ($HMM1, $HMM2) = $self->hmms;
    $out .= sprintf ("%7s %7d %7d %7d         1 %7s %7d %7d %7d %7s %7s %7s\n", $HMM1->accession, $self->{'start1'}, $self->{'end1'}, $self->{'length1'}, $HMM2->accession, $self->{'start2'}, $self->{'end2'}, $self->{'length2'}, $self->score, $self->{'simple'}, $self->{'reverse'});
    $out .= "\n".$self->{'string1'}."\n".$self->{'string2'}."\n";
    $out .= "\n# END\n";
}

sub algorithm 
{
    my $self = shift;
    return 'forward/backward' if $self->{'algorithm'} && $self->{'algorithm'} eq 'forward';
    return $self->{'algorithm'};
}

sub version
{
    my $self = shift;
    return $self->{'version'};
}

sub mode
{
    my $self = shift;
    return $self->{'mode'};
}

sub hmms
{
    my $self = shift;
    return ($self->{'HMM1'}, $self->{'HMM2'});
}

sub length
{
    my $self = shift;
    my $alignment = $self->{'string1'};
    $alignment =~ s/~//gi if $alignment;
    return length($alignment);
}

sub score
{
    my $self = shift;
    return $self->{'co-emis'} if $self->{'co-emis'};
    return "NA";
}

# A private function to populate the member variables in this object with values from a prc file.
sub _parseFile
{
    my ($self, $file) = @_;
    
    ($self->{'version'}) = $file =~ /PRC ([\d\.]+)/;
    ($self->{'command'}) = $file =~ /Command\s+:\s+(.+)\s*\n/;
    ($self->{'algorithm'}) = $file =~ /Algorithm\s+:\s+(\w+)/;
    ($self->{'mode'}) = $file =~ /Align mode\s+:\s+(\S+)/;
    ($self->{'date'}) = $file =~ /Started\s+:\s+(.+)\s*\n/;
    ($self->{'stop'}) = $file =~ /Simple stop\s+:\s+(\S+)/;
    
    #A very large pattern that breaks the prc-style output (as opposed to "sam1/2" styles) of prc into its components
    $file =~ /#\s+(hmm.+)\s*\n(.+)\n/i;
    @{$self}{split(/\s+/, "\L$1")} = split(/\s+/, $2); # Split the values into a hash with corresponding names from the line above
    
    #The alignment strings
    $file =~ /\n\s*\n([MDI~]*)\n([MDI~]*)\n/;
    $self->{'string1'} = $1;
    $self->{'string2'} = $2;
    # Remove trailing/heading gaps
    if($self->{'string1'} =~ /^(D+)/)
    {
        $self->{'start1'} += CORE::length($1);
        $self->{'string1'} =~ s/^[D~]+//;
        $self->{'start2'}++;
    }
    $self->{'end1'} -= CORE::length($1) if $self->{'string1'} =~ /(D+)$/;
    $self->{'string1'} =~ s/[D~]+$//;
    
    if($self->{'string2'} =~ /^(D+)/)
    {
        $self->{'start2'} += CORE::length($1);
        $self->{'string2'} =~ s/^[D~]+//;
        $self->{'start1'}++;
    }
    $self->{'end2'} -= CORE::length($1) if $self->{'string2'} =~ /(D+)$/;
    $self->{'string2'} =~ s/[D~]+$//;
}

# This is a private method that returns the raw logo part of the image
sub _createLogo 
{
    my ($self, $ICM, $HPM, $width, $alphabet, $startPos, $xSize, $ySize, $alString, $greyscale, $bottomNumbers) = @_;
    
    $width = pdl $width; # Work on Copy!
    
    #Imager-Objects
    my $black=NC(0,0,0);
    my $transBlack=Imager::Color->new(r=>0, g=>0, b=>0, a=>128);
    my $white=NC(255,255,255);
    
    #fill the bars for insert states with either a solid color or slashes
    my $hpFill = $greyscale?(Imager::Fill->new(hatch=>'dots16', fg=>$black, bg=>$white)):(Imager::Fill->new(solid => NC('#FF6666')));
    my $contribFill=$greyscale?(Imager::Fill->new(hatch=>'dots4', fg=>$black, bg=>$white)):(Imager::Fill->new(solid => NC('#FFCCCC')));
    my $unalignedFill = Imager::Fill->new(solid=>$transBlack, combine=>'normal');

    my @alphabet = @{$alphabet};
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
                             'G' => NC('#FFFF00'),
                             'T' => NC('#FF0000'),
                             'U' => NC('#FF0000')
                             );
        }
    }
    my $normFont=NF(file=>$REGFONT, type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";
    my $titleFont=NF(file=>$BOLDFONT, type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";
    my $vertFont=NF(file=>$REGFONT, type=>'ft2') or die "Couldn't create font: $Imager::ERRSTR\n";

    my $normFontSize = 15; # Should this be dynamical?
    my $smallFontSize = 12;
    my $tinyFontSize = 10;

    #make new image, fill white
    my $i=Imager->new(xsize=>$xSize, ysize=>$ySize, channels=>($greyscale?1:4)); # destination image
    $i->box(filled=>1,color=>$white);                     # fill with background color

    #initialize the most important variables first
    my $maxInfContent = _max(list $ICM->xchg(0,1)->sumover()); # How high should the y-axis be?
    my $motifSize = $ICM->getdim(0); # And how many states do we have?

    #get the original HPM-Values for Insert states
    my $HPV = ones($motifSize);
    $HPV->slice('1:'.($motifSize-1).':2') .= $HPM->slice(':,1');

    my $bottomMargin = $ySize;
    $bottomMargin -= $bottomNumbers?15:5;
    my $topMargin = ($bottomNumbers?5:15);
    my $leftMargin = 15;
    my $rightMargin = $xSize - 5;
    
    my $logoHeight = $ySize - 20; # This should be dynamic, too. 15+5 = 20, thats the sum of the upper+lower reserved area for numbers etc.
    my $xScaleFactor = ($rightMargin - $leftMargin)/$width->sumover(); # Scale to fit the longer HMM

    #scale width to target-width 
    $width *= $xScaleFactor;
    
    #scale HP to target-width 
    $HPV *= $xScaleFactor;    
    
    # Draw Y-Axes
    $i->line(color=>$black, 
             x1=> $leftMargin,
             x2=> $leftMargin,
             y1=> $bottomMargin,
             y2=> $topMargin);

    # Draw Y-Axis ticks
    my $yStep = $logoHeight / PDL::Math::ceil($maxInfContent);
    for(my $k=0; $k<=PDL::Math::ceil($maxInfContent); ++$k){
        $i->line(color=>$black,
                x1=>$leftMargin-3,
                x2=>$leftMargin+3,
                y1=>$bottomMargin - $k*$yStep,
                y2=>$bottomMargin - $k*$yStep);
    
        $i->string(font=>$normFont,
                   string=>"$k ",
                   x=>$leftMargin - 11, 
                   y=>$bottomMargin - $k*$yStep + 4,
                   size=>$tinyFontSize,
                   color=>$black,
                   aa=>1);
    }

    # Draw the X-axes
    $i->line(color=>$black, 
             x1=> $leftMargin,
             x2=> $rightMargin,
             y1=> $bottomMargin,
             y2=> $bottomMargin);
    
    # Set up a boolean vector that determines whether a state is to be greyed out or not
    my $isAligned = zeroes($motifSize);
    
    # Remove the gaps from the string. They don't alter anything concerning this HMM
    $alString =~ s/~//g;
    
    #The first letter of the alignment string will almost certainly be a match. Then we don't need to do anything, as $isAligned->at(0) == 0 already!
    my $curState = 0;
    #  The other possibility is that it's an insert, for which we take care here:
    $isAligned->slice('0') .= 1 if ( unpack('A1', $alString) eq 'I' );
    # If we have a leading Delete. Don't tell me this is bollocks, cause I know that it is, but Martin likes the idea
    if ( unpack('A1', $alString) eq 'D' )
    {
        $isAligned->slice('0:1') .= 1;
        $curState += 1;
    }
    for (my $curChar=1; $curChar < CORE::length($alString); ++$curChar)
    {
        my ($lastChar, $alignedChar) = unpack('x'.($curChar-1).' A1 A1', $alString);
        
        # Now, we black out the unmatched fields
        if ($alignedChar eq 'M')
        {
            if ($lastChar eq 'M')
            {
                $curState += 2; # We progres to the next match state, skipping one insert!
                # We have to grey out the last insert state
                $isAligned->slice($curState-1) .= 1;
            }
            elsif ($lastChar eq 'I')
            {
                $curState += 1; # We progres to the next match state
            }
            elsif ($lastChar eq 'D')
            {
                $curState += 1; # We progres to the next match state
            }
        }
        elsif ($alignedChar eq 'I')
        {
            if ($lastChar eq 'M')
            {
                $curState += 1; # We proceed into this Insert state
            }
        }
        elsif ($alignedChar eq 'D')
        {
            if ($lastChar eq 'M')
            {
                $curState += 3;
                $isAligned->slice(($curState - 2).':'.$curState) .= 1; # The next 3 states (IMI) are skipped!
            }
            elsif ($lastChar eq 'D')
            {
                $curState += 2;
                $isAligned->slice(($curState - 1).':'.$curState) .= 1; # The next 2 states (MI) are skipped!
            }
        }
        else
        {
            # This should never happen!
            warn("FATAL: Received invalid alignment letter: ", $alignedChar);
            return 0;
        }
    }

    # Grey out the trailing insert if we end in a match!
    ($isAligned->slice(++$curState) .= 1) if (unpack('x'.(CORE::length($alString)-1).' A1', $alString) eq 'M');
    
    # Draw the Logo
    # Boolean switch that indicates whether we're in a match or an insert state...
    my $isInsert = 0;
    for(my $k=0; $k<$motifSize; ++$k)
    {
        # Reset the infomation content hash
        my %infContent = ();
        #get all the values for one column from the matrix
        my @infValues = list $ICM->slice($k);

        #fill infContent-Hash with Values
        my $l=0;
        foreach my $key (@alphabet){
            unless(($infValues[$l]*$yStep)<1){ # if the character would print lower than 1 pixel, discard it
                $infContent{$key} = sprintf('%.16f', $infValues[$l]);
            }
            ++$l; # this is the index of the information-content of the next character in the infValues vector
        }
        
        #define drawing order by descending infContent
        my @drawOrder = sort {$infContent{$a}<=>$infContent{$b}} keys(%infContent); # Sort the letters to be printed by decreasing inf-content
    
        #reset the height to the bottom Margin
        my $sumHeight = $bottomMargin;
        
        if($isInsert)
        {
            # We're in an insert, so let's print the state-number, draw the x-tick and create the red/pink blocks
            #draw tick and number
            $i->line(color=>$black,
                    x1  => $leftMargin + $width->slice("0:$k")->sumover(),
                    x2  => $leftMargin + $width->slice("0:$k")->sumover(),
                    y1  => ($bottomNumbers?$bottomMargin:$topMargin) - 3,
                    y2  => ($bottomNumbers?$bottomMargin:$topMargin) + 3);
                    
            $i->string(font => $normFont,
                   string   => (($startPos + $k+1)/2).' ',
                   x        => $leftMargin + $width->slice("0:$k")->sumover() - 
                               0.5 * $width->slice(($k-1).':'.$k)->sumover() - 
                               0.5 * $normFont->bounding_box(string=>($startPos + $k+1)/2)->total_width, 
                   y        => $bottomNumbers?($bottomMargin + 13):($topMargin-3),
                   size     => $tinyFontSize,
                   color    => $black,
                   aa       =>  1);
                   
            #draw bg for contribution
            $i->box(fill=>$contribFill, #use the fill selected earlier
                    xmin=>$leftMargin + $width->slice("0:".($k-1))->sumover(),
                    ymin=>$topMargin,
                    xmax=>$leftMargin + $width->slice("0:$k")->sumover(),
                    ymax=>$bottomMargin-1);
                    
            #draw bg for hitting probability
            $i->box(fill=>$hpFill, #use the fill selected earlier
                    xmin=>$leftMargin + $width->slice("0:".($k-1))->sumover(),
                    ymin=>$topMargin,
                    xmax=>$leftMargin + $width->slice("0:".($k-1))->sumover() + $HPV->at($k),
                    ymax=>$bottomMargin-1);
            $isInsert = 0;
        }
        else
        { # We're in a match, don't do anything special...
            $isInsert = 1;
        }
        
        foreach my $key (@drawOrder){
#            #get the size of the character to draw
            my $bbox = $normFont->bounding_box(string=>$key, size=>$infContent{$key}*$yStep, sizew=>$width->at($k,0));
#            #there is a difference between target height and real drawn height, because many characters do not use all the space provided.
#            #To solve that issue, we calculate a factor to multiply the target height with, so it compensates the difference in sizes.
            my $scaleFactor = $infContent{$key}* $yStep / $bbox->text_height * 0.98;
            
            # We can't use total_width because it has a nasty bug! Instead, we use adv_width - neg_width
            my $scaleFactorW = $width->at($k,0)/($bbox->advance_width - $bbox->neg_width);
            
            if($infContent{$key}){
                # Plot the next letterin the stack
                $i->string(font=>$normFont,
                           string=>$key,
                            # We want centered text; thus, we start from the middle of the column and go left by half of the positiv width of the letter.
                           x=>$leftMargin + ($k?($width->slice("0:".($k-1))->sumover()):0) + (0.5*$width->at($k,0)-0.5*$bbox->advance_width*$scaleFactorW),
                           y=>$sumHeight + $bbox->descent*$scaleFactor,
                           size=>$infContent{$key} * $yStep * $scaleFactor,
                           sizew=>$width->at($k,0)*$scaleFactorW*0.98,
                           color=>$letterColors{$key},
                           aa=>1);
                $sumHeight -= $bbox->text_height * $scaleFactor; #$infContent{$key}*$yStep - ($bbox->descent*$scaleFactor); # Old style...
            }
        }
        
        # Draw black rect around the unaligned states
        if ($isAligned->at($k))
        {
            $i->box(fill=>$unalignedFill, # This is supposed to be 50% transparent
                    xmin=>$leftMargin + ($k?($width->slice("0:".($k-1))->sumover()):0),
                    ymin=>$topMargin,
                    xmax=>$leftMargin + $width->slice("0:$k")->sumover(),
                    ymax=>$bottomMargin-1);
        }
        else
        {
            $i->box(color=>$black, #fill=>$unalignedFill, # This is supposed to be 50% transparent
                    xmin=>$leftMargin + ($k?($width->slice("0:".($k-1))->sumover()):0),
                    ymin=>$topMargin,
                    xmax=>$leftMargin + $width->slice("0:$k")->sumover() - 1,
                    ymax=>$bottomMargin);
        }
    }
    
    return $i;
}

# Return the maximum value in a list
sub _max {
    my @list = @_;
    warn "usage: max @list where @list !>= 2" unless @list >= 2;

    my $max = shift @list;
    foreach my $elem (@list)
    {
        $max = $elem if $elem > $max;
    }
    return $max;
}

return 1; 

=head1 BUGS

This package is still in experimental stage; if you witness malfunction please contact the author.
The author will not take any responsibility for loss of time, data, money or other occuring problems.
Use it at your own risk!

=head1 AUTHOR

Benjamin Schuster-Boeckler E<lt> bsb(at)sanger.ac.uk E<gt>

=head1 SEE ALSO

L<HMM::Profile> L<HMM::Utilities>, L<PDL>

=cut
