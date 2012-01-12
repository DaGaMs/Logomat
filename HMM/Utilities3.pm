package HMM::Utilities3;
use vars qw(@ISA @EXPORTER @EXPORT_OK);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = ();

use PDL::LiteF;

# Utilities3.pm
# subroutines necessary to convert a HMMER3 *.hmm file
# into its probability format.
# modified to use PDL

my $MININTVAL = -987654321;
my $INTSCALE = 1000;

sub Prob2Log {    #Convert probabilities to scores
    my $prob = shift;
    my $mask = shift;
    return -1 * $prob->log()->badmask($mask);
}
sub Log2Prob {    #Convert scores (a piddle, normally) to probability
    my $score = shift;
    return exp(-1 * $score);
}

sub Ascii2Prob {    #Convert input values (* & nums) to Probs
    my $ascii = shift;
    chomp($ascii);
    $ascii =~ s/\*/$MININTVAL/g;
    my $score = pdl split(/\s+/, $ascii); #return a piddle
    return &Log2Prob($score);
}

sub Prob2Ascii {    #Converts probs-vectors to single line, space separated strings
    return Log2Ascii(Prob2Log(@_));
}

sub Ascii2Log {    #Convert input values (* & nums) to Scores
    my $ascii = shift;
    chomp($ascii);
    $ascii =~ s/\*/$MININTVAL/g;
    return pdl split(/\s+/, $ascii); #return a piddle
}

sub Log2Ascii {    #Converts Score-vectors to single line, space separated strings
    my $score = shift;
    my $out = '';
    foreach my $number (list($score))
    {
        if ($number <= $MININTVAL)
        {
        $out .= sprintf(' %7s', '*'); # Create a string where numbers are right padded with varying distance
        }
        else
        {
        $out .= sprintf(' %01.5f', $number); # Create a string where numbers are right padded with varying distance
        }
        
    }
    return $out;
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

=item Prob2Log($prob, $nullModel)

 Usage     : my $score = Prob2Log($prob, $nullModel)
 Function: convert probability-piddles to score-piddles
 Returns : a new piddle

=item Log2Prob($score, $nullModel)

 Usage     : my $prob = Prob2Log($score, $nullModel)
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

=item Ascii2Log($ascii)

 Usage     : my $score = Ascii2Log($ascii)
 Function: convert perl scalars to score-piddles
 Returns : a new piddle

=item Log2Ascii($score)

 Usage     : my $ascii = Ascii2Prob($score)
 Function: convert score-piddles to string of integers
 Returns : a string

=head1 BUGS

This package is still in experimental stage; if you witness malfunction please contact the author.
The author will not take any responsibility for loss of time, data, money or other occuring problems.
Use it at your own risk!

=head1 AUTHOR

Benjamin Schuster-Boeckler E<lt> boeckler(at)molgen.mpg.de E<gt>

=head1 SEE ALSO

L<HMM::Profile>, L<HMM>, L<PDL>

=cut
