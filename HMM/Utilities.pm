package HMM::Utilities;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(toICM toHPM);

use PDL::LiteF;



# Utilities.pm
# subroutines necessary to convert an *.hmm file
# into its probability format.
# modified to use PDL

# was: hmmfunc.pl from hmmerviewer

my $LBASE2 = 1.44269504;
my $EBASE2 = 0.69314718;
my $INTSCALE = 1000;

sub Prob2Score {    #Convert probabilities to scores
    my $prob = shift;
    my $nul_model = shift;
    my $mask = shift;
    return (0.5 + $INTSCALE*&sreLOG2($prob/$nul_model, $mask))->floor();
}
sub Score2Prob {    #Convert scores (a piddle, normally) to probability
    my $score = shift;
    my $nul_model = shift;
    return ($nul_model * 2**($score/$INTSCALE));
}

sub Ascii2Prob {    #Convert input values (* & nums) to Probs
    my $ascii = shift;
    my $nul_model = shift;
    chomp($ascii);
    $ascii =~ s/\*/-987654321/g;
    my $score = transpose(pdl split(/\s+/, $ascii)); #return a piddle
    return &Score2Prob($score, $nul_model);
}

sub Prob2Ascii {    #Converts probs-vectors to single line, space separated strings
    return Score2Ascii(Prob2Score(@_));
}

sub Ascii2Score {    #Convert input values (* & nums) to Scores
    my $ascii = shift;
    chomp($ascii);
    $ascii =~ s/\*/-987654321/g;
    return transpose(pdl split(/\s+/, $ascii)); #return a piddle
}

sub Score2Ascii {    #Converts Score-vectors to single line, space separated strings
    my $score = shift;
    my $out = '';
    foreach my $number (list($score))
    {
        $out .= sprintf(' %6d', $number); # Create a string where numbers are right padded with varying distance
    }
    $out =~ s/-987654321/     */g;
    return $out;
}

sub sreLOG2 {        #Log base 2 as used in HMMer
    my $x = shift;
    my $mask = shift || 0;
    return (log($x)*$LBASE2)->badmask($mask);
}

#create an information content matrix from the emission-probabilities
sub toICM{
    my $prob = shift;
    my $nul_model = shift;
    my $height_logodds = shift; # boolean.  If 0, then use emission probability for height (which was the default) ,  TJW

    my $xsize = $prob->getdim(0);
    my $ysize = $prob->getdim(1);
    my $flat = zeroes($xsize*2, $ysize);

    $flat->slice('0:'.($xsize*2-2).':2,:') .= $prob->slice(':,:,0');
    $flat->slice('1:'.($xsize*2-1).':2,:') .= $prob->slice(':,:,1');

    #get back to scores
    my $ICM = $flat * &sreLOG2($flat/$nul_model);

    #normalize to column-height
    if ($height_logodds) {
         #new style, which splits the information content height according to log odds
         my $odds =  &sreLOG2($flat/$nul_model);
         $odds    =  $odds->setbadif( $odds <= 0);  # only include entries with positive log-odds
         $odds    =  $odds / $odds->xchg(0, 1)->sumover();
         $ICM     =  $odds * $ICM->xchg(0, 1)->sumover();
    } else {
         #original style, which splits the information content height according to emission probability
         $ICM     = $flat * $ICM->xchg(0, 1)->sumover();
    }

    return $ICM;
}


#calculate the probabilities of entering a state (results in column-width of the logo)
sub toHPM{
    my $startTransitions = shift;
    my $transitions = shift;

    my $xSize = $transitions->getdim(0);
    my $width = zeroes($xSize, 3);
    if ($startTransitions->getdim(1) == 3)
    {
        # HMMER2
        $width->slice('0,:') .= $startTransitions; #initialize with start-probabilities
        $width->slice('0,1') .= $startTransitions->at(0,0) * $transitions->at(0, 1); #initialize first insert
    }
    elsif ($startTransitions->getdim(1) == 7)
    {
        # HMMER3
	# For the first match state, there is the possibility of going through insert state 0, too
	# In HMMER2, there was no insert state 0
	# This means that currently, Insert state 0 will not be displayed!
	$width->slice('0,0') .= $startTransitions->at(0, 1) + $startTransitions->at(0, 0);
	$width->slice('0,1') .= $width->slice('0,0') * $transitions->at(0, 1);
	$width->slice('0,2') .= $startTransitions->at(0, 2);
    }
    for(my $i=1; $i<$xSize; ++$i){
        #Match
        $width->slice("$i,0") .= $width->at($i-1, 1) + $width->at($i-1, 0)*$transitions->at($i-1, 0) + $width->at($i-1, 2)*$transitions->at($i-1, 5);
        #Insert
        $width->slice("$i,1") .= ($width->at($i, 0) * $transitions->at($i, 1));
        #Delete
        $width->slice("$i,2") .= $width->at($i-1, 0)*$transitions->at($i-1, 2) + $width->at($i-1, 2)*$transitions->at($i-1, 6);
    }
    return $width;
}

=head1 NAME

HMM::Utilities - Utility functions for HMMs

=head1 SYNOPSIS

use HMM::Utilities;

=head1 DESCRIPTION

This module contains helpfull methods for converting scores to probabilities and vice versa
as well as creation of information content and hitting probability matrices.

=head1 METHODS

=over 4

=item Prob2Score($prob, $nullModel)

 Usage     : my $score = Prob2Score($prob, $nullModel)
 Function: convert probability-piddles to score-piddles
 Returns : a new piddle

=item Score2Prob($score, $nullModel)

 Usage     : my $prob = Prob2Score($score, $nullModel)
 Function: convert score-piddles to probability-piddles
 Returns : a new piddle

=item Ascii2Prob($ascii, $nullModel)

 Usage     : my $prob = Ascii2Prob($ascii, $nullModel)
 Function: convert perl scalars to probability piddles
 Returns : a new piddle

=item Prob2Ascii($prob, $nullModel)

 Usage     : my $ascii = Prob2Ascii($prob)
 Function: convert probability piddles to string of floats
 Returns : a string

=item Ascii2Score($ascii)

 Usage     : my $score = Ascii2Score($ascii)
 Function: convert perl scalars to score-piddles
 Returns : a new piddle

=item Score2Ascii($score)

 Usage     : my $ascii = Ascii2Prob($score)
 Function: convert score-piddles to string of integers
 Returns : a string

=item sreLog2($val)

 Usage     : my $odds = sreLog2($piddle)
 Function: do Log to base 2 the hmmer way on a piddle
 Returns : a new piddle

=item toICM($prob, $nullModel)

 Usage     : my $ICM = toICM($prob, $nullModel)
 Function: generate a new information content matrix from a probability matrix
 Returns : a new 2d-piddle

=item toHPM($prob, $nullModel)

 Usage     : my $HPM = toHPM($startTransitions, $transitions)
 Function: generate a new hitting-probability matrix
 Returns : a new piddle

 =back

=head1 BUGS

This package is still in experimental stage; if you witness malfunction please contact the author.
The author will not take any responsibility for loss of time, data, money or other occuring problems.
Use it at your own risk!

=head1 AUTHOR

Benjamin Schuster-Boeckler E<lt> boeckler(at)molgen.mpg.de E<gt>

=head1 SEE ALSO

L<HMM::Profile>, L<HMM>, L<PDL>

=cut
