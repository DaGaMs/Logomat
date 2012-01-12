=head1 NAME

HMM - representation of hidden markov moddels for computational biology

=head1 SYNOPSIS

  use HMM;
    
  my $hmm = HMM->new(-hmmerfile=>$file) || die("Couldn't open $file!\n"); #try to create a new hmm-object from the file
    
  print $hmm->name(); #output the name of the model
  print $hmm->transitions(); #show the matrix of transition probabilities

=head1 DESCRIPTION

This module implements several basic access methods to HMMER models. On creation, it parses a HMMER file either
provided as a filename, a filehandle or as a PFAM id. The user can then retrieve the information contained in the
file through the methods described below. HMM relies strongly on the Perl Data Language, which makes it possible
to manage matrices in a MATLAB-like manner. If you are not familiar with PDL read the documentation provided with
it: L<PDL>.

=head1 METHODS

=cut

package HMM;

use strict;
use PDL::LiteF;

our $VERSION = '0.81';

=head2 new(%args)

 Usage      : my $pHMM = HMM->new(%args)
 Function   : constructor for the HMM object
 Returns    : a new HMM object
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
    bless($self, ref($class) || $class);
    $self->{'alphabet'} = $args{-alphabet} if ($args{-alphabet});
    $self->{'alphabet'} = [qw(A C D E F G H I K L M N P Q R S T V W Y)] unless $self->{'alphabet'};
    $self->{'transitions'} = $args{-transitions} if (ref $args{-transitions});
    $self->{'startTransitions'} = $args{-startTransitions} if (ref $args{-startTransitions});
    $self->{'startTransitions'} = transpose(pdl( qw(1 0 0) )) unless ref $self->{'startTransitions'};
    $self->{'emissions'} = $args{-emissions} if (ref $args{-emissions});
    $self->{'nullEmissions'} = $args{-nullEmissions} if (ref $args{-nullEmissions});
    my $prob = 1 / scalar(@{$self->{'alphabet'}});
    unless (ref $self->{'nullEmissions'})
    {
        my @unif = split ' ', "$prob " x scalar(@{$self->{'alphabet'}});
        $self->{'nullEmissions'} = transpose(pdl(@unif));
    }
    $self->{'accession'} = $args{-accession} if ($args{-accession});
    $self->{'version'} = $args{-version} || $VERSION;
    $self->{'seqNumber'} = $args{-seqNumber} if ($args{-seqNumber});
    $self->{'modelName'} = $args{-name} if ($args{-name});
    $self->{'description'} = $args{-description} if ($args{-description});

    return $self;
}

=head2 version()

 Usage      : my $text = $pHMM->version()
 Returns    : the version-number of the HMMer implementation

=cut

sub version {
    my $self = shift;
    return $self->{'version'};
}

=head2 name()

 Usage      : my $text = $pHMM->name()
 Returns    : the name of the represented family

=cut

sub name : lvalue {
    my $self = shift;
    $self->{'modelName'};
}

    
=head2 alphabet()

 Usage      : my @alphabet = $pHMM->alphabet()
 Returns    : either [A,C,G,T] or the one-letter codes of the 20 AAs

=cut

sub alphabet : lvalue {
    my $self = shift;
    $self->{'alphabet'};
}

    
=head2 length()

 Usage      : my $length = $pHMM->length()
 Returns    : the number of columns the alignment consists of

=cut

sub length {
    my $self = shift;
    return $self->{'length'};
}

=head2 seqNumber()

 Usage      : my $number = $pHMM->seqNumber()
 Returns    : the number of sequences the profile was built with

=cut

sub seqNumber {
    my $self = shift;
    return $self->{'seqNumber'};
}

=head2 specialTransitions()

 Usage      : my $piddle = $pHMM->specialTransitions()
 Returns    : the transition probabilities for the N,E,C and J states as defined by the HMMer-Userguide

=cut

sub specialTransitions : lvalue {
    my $self = shift;
    $self->{'specialTransitions'};
}

=head2 nullTransitions()

 Usage      : my $piddle = $pHMM->nullTransitions()
 Returns    : two transition probabilities for staying in or leaving the background state

=cut

sub nullTransitions : lvalue {
    my $self = shift;
    $self->{'nullTransitions'};
}

=head2 nullEmissions()

 Usage      : my $piddle = $pHMM->nullEmissions()
 Returns    : background emission probabilities for the used alphabet

=cut

sub nullEmissions : lvalue {
    my $self = shift;
    $self->{'nullEmissions'};
}

=head2 startTransitions()

 Usage      : my $piddle = $pHMM->startTransitions()
 Returns    : transition probabilities for B->M1 (B->I) B->D1

=cut

sub startTransitions : lvalue {
    my $self = shift;
    $self->{'startTransitions'};
}

=head2 transitions()

 Usage      : my $piddle = $pHMM->transitions()
 Returns    : transition probabilities for every column of the profile; Transitions are:
     M->M  M->I  M->D  I->M  I->I  D->M  D->D  B->M  M->E
    The dimensions are $pHMM->length() x 9

=cut

sub transitions : lvalue {
    my $self = shift;
    $self->{'transitions'};
}

=head2 emissions()

 Usage      : my $piddle = $pHMM->emissions()
 Returns    : 2-dimensional piddle of emission probabilities for each M and I state;
     $pHMM->emissions->slice(':,:,0') gives M, ...slice(':,:,1') gives I

=cut

sub emissions : lvalue {
    my $self = shift;
    $self->{'emissions'};
}

=head2 specialTransitionScores()

 Usage      : my $piddle = $pHMM->specialTransitions()
 Returns    : the transition scores for the N,E,C and J states as defined by the HMMer-Userguide

=cut

sub specialTransitionScores : lvalue {
    my $self = shift;
    $self->{'specialTransitionScores'};
}

=head2 nullTransitionScores()

 Usage      : my $piddle = $pHMM->nullTransitions()
 Returns    : two transition scores for staying in or leaving the background state

=cut

sub nullTransitionScores : lvalue {
    my $self = shift;
    $self->{'nullTransitionScores'};
}

=head2 nullEmissionScores()

 Usage      : my $piddle = $pHMM->nullEmissions()
 Returns    : background emission scores for the used alphabet

=cut

sub nullEmissionScores : lvalue {
    my $self = shift;
    $self->{'nullEmissionScores'};
}

=head2 startTransitionScores()

 Usage      : my $piddle = $pHMM->startTransitions()
 Returns    : transition scores for B->M1 (B->I) B->D1

=cut

sub startTransitionScores : lvalue {
    my $self = shift;
    $self->{'startTransitionScores'};
}

=head2 transitionScores()

 Usage      : my $piddle = $pHMM->transitionScores()
 Returns    : transition scores for every column of the profile; Transitions are:
    M->M  M->I  M->D  I->M  I->I  D->M  D->D  B->M  M->E
    The dimensions are $pHMM->length() x 9

=cut

sub transitionScores : lvalue {
    my $self = shift;
    $self->{'transitionScores'};
}

=head2 emissionScores()

 Usage      : my $piddle = $pHMM->emissionScores()
 Returns    : 2-dimensional piddle of emission scores for each M and I state;
     $pHMM->emissions->slice(':,:,0') gives M, ...slice(':,:,1') gives I

=cut

sub emissionScores : lvalue {
    my $self = shift;
    $self->{'emissionScores'};
}

=head2 accession()

 Usage      : my $acc = $pHMM->accession()
 Returns    : accession as string

=cut

sub accession : lvalue {
    my $self = shift;
    $self->{'accession'};
}

=head2 description()

 Usage      : my $desc = $pHMM->accession()
 Returns    : description as string

=cut

sub description : lvalue {
    my $self = shift;
    $self->{'description'};
}

=head2 evidence()

 Usage      : my $evd = $pHMM->evidence()
 Returns    : lambda and nu parameters of E-Value distribution of this HMM, if calibrated

=cut

sub evidence : lvalue {
    my $self = shift;
    $self->{'evidence'};
}

=head2 state2column($state)

 Usage      : my $col = $pHMM->state2column($state)
 Returns    : column in the seed alignment that corresponds to the given state index (index starting with 1)

=cut

sub state2column {
    my $self = shift;
    my $state = shift;
    return undef unless $state;
    return $self->{'state2column'}->[$state-1];
}

=head2 column2state($column)

 Usage      : my $state = $pHMM->column2state($column)
 Returns    : state in the model that corresponds to the column in the seed alignment

=cut

sub column2state {
    my $self = shift;
    my $col = shift;
    return undef unless $col;
    return undef if $col > $self->{'state2column'}->[@{$self->{'state2column'}} - 1];
    my $state;
    map {++$state if $col >= $_;} @{$self->{'state2column'}};
    return $state;
}

1;

=head1 BUGS

This package is still in experimental stage; if you witness malfunction please contact the author.
The author will not take any responsibility for loss of time, data, money or other occuring problems.
Use it at your own risk!

=head1 AUTHOR

Benjamin Schuster-Boeckler E<lt> bsb(at)sanger.ac.uk E<gt>

=head1 SEE ALSO

L<HMM::Profile> L<HMM::Utilities>, L<PDL>

=cut
