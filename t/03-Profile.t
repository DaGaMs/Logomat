#test-script for HMM::Profile
use Test::More tests => 72;

BEGIN { use_ok(HMM::Profile); }

#
# Test if the hmmer2 parser works correctly
#
my $file = 't/Piwi.hmm';

my $hmm = HMM::Profile->new(-hmmerfile => $file);

#Test object creation
ok(defined $hmm, "HMM::Profile->new returned value");
isa_ok($hmm, "HMM::Profile");

# Test annotation
is($hmm->name, 'Piwi', "Name correctly parsed");
is($hmm->version, '2.0', 'Version correctly parsed');
is($hmm->accession, 'PF02171.9', "Accession correctly parsed");
is($hmm->description, 'Piwi domain', "Description correctly parsed");
cmp_ok($hmm->length, '==', 346, 'Length correctly parsed');
cmp_ok($hmm->seqNumber, '==', 21, 'Number of seqs correctly parsed');
is_deeply($hmm->alphabet, [qw(A C D E F G H I K L M N P Q R S T V W Y)], "Alphabet correctly parsed");

# Test model
cmp_ok($hmm->nullEmissionScores->getdim(1), '==', 20, 'Null emission scores of right size');
cmp_ok($hmm->nullEmissions->getdim(1), '==', 20, 'Null emissions of right size');
cmp_ok($hmm->nullTransitionScores->getdim(1), '==', 2, 'Null transition scores of right size');
cmp_ok($hmm->nullTransitions->getdim(1), '==', 2, 'Null transitions of right size');
cmp_ok($hmm->emissions->getdim(0), '==', 346, 'Emissions have right number of states');
cmp_ok($hmm->emissions->getdim(1), '==', 20, 'Emissions have right number of residues');
cmp_ok($hmm->emissions->getdim(2), '==', 2, 'Emissions have both match and insert dimensions');
cmp_ok($hmm->transitions->getdim(0), '==', 346, 'Transitions have right number of states');
cmp_ok($hmm->transitions->getdim(1), '==', 9, 'Transitions have right number of transitions');

cmp_ok($hmm->nullEmissionScores->at(0,0), '==', 595, 'Right null emission score for A');
cmp_ok($hmm->nullEmissionScores->at(0,1), '==', -1558, 'Right null emission score for C');
cmp_ok($hmm->nullEmissionScores->at(0,18), '==', -1998, 'Right null emission score for W');
cmp_ok($hmm->nullEmissionScores->at(0,19), '==', -644, 'Right null emission score for Y');
is(sprintf('%1.5f', $hmm->nullEmissions->at(0,0)), '0.07552', 'Right null emission prob for A');
is(sprintf('%1.5f', $hmm->nullEmissions->at(0,4)), '0.04078', 'Right null emission prob for F');
is(sprintf('%1.2f', $hmm->nullEmissions->squeeze()->sumover()), '1.00', 'Cross sum of null emissions is 1');

cmp_ok($hmm->emissionScores->at(0, 0, 0), '==', -1981, 'Right match emission score for A at 1');
cmp_ok($hmm->emissionScores->at(10, 1, 0), '==', -3254, 'Right match emission score for C at 11');
cmp_ok($hmm->emissionScores->at(300, 15, 0), '==', 15, 'Right match emission score for S at 301');
cmp_ok($hmm->emissionScores->at(345, 19, 0), '==', -2775, 'Right match emission score for Y at 346');

cmp_ok($hmm->emissionScores->at(0, 0, 1), '==', -149, 'Right insert emission score for A at 1');
cmp_ok($hmm->emissionScores->at(10, 1, 1), '==', -500, 'Right insert emission score for C at 11');
cmp_ok($hmm->emissionScores->at(300, 15, 1), '==', 359, 'Right insert emission score for S at 301');
cmp_ok($hmm->emissionScores->at(345, 19, 1), '==', -987654321, 'Right insert emission score for Y at 346');

is(sprintf('%1.5f', $hmm->emissions->at(0, 0, 0)), '0.01913', 'Right match emission prob for A at 1');
is(sprintf('%1.5f', $hmm->emissions->at(14, 4, 0)), '0.00778', 'Right match emission prob for F at 15');
is(sprintf('%1.2f', $hmm->emissions->slice('250,:,0')->squeeze()->sumover()), '1.00', 'Cross sum of match emissions at 251 is 1');
is(sprintf('%1.2f', $hmm->emissions->slice('133,:,1')->squeeze()->sumover()), '1.00', 'Cross sum of insert emissions at 134 is 1');

############### Now test HMMER3 parsing ###############

$file = 't/Pkinase.h3m';

$hmm = HMM::Profile->new(-hmmerfile => $file);

#Test object creation
ok(defined $hmm, "HMM::Profile->new returned value");
isa_ok($hmm, "HMM::Profile");

# Test annotation
is($hmm->name, 'Pkinase', "Name correctly parsed");
is($hmm->version, '3/a', 'Version correctly parsed');
is($hmm->accession, 'PF00069.16', "Accession correctly parsed");
is($hmm->description, 'Protein kinase domain', "Description correctly parsed");
cmp_ok($hmm->length, '==', 260, 'Length correctly parsed');
cmp_ok($hmm->seqNumber, '==', 54, 'Number of seqs correctly parsed');
is_deeply($hmm->alphabet, [qw(A C D E F G H I K L M N P Q R S T V W Y)], "Alphabet correctly parsed");

# Test model
cmp_ok($hmm->nullEmissionScores->getdim(1), '==', 20, 'Null emission scores of right size');
cmp_ok($hmm->nullEmissions->getdim(1), '==', 20, 'Null emissions of right size');
# Null transitions don't exist any more in HMMER3
#cmp_ok($hmm->nullTransitionScores->getdim(1), '==', 2, 'Null transition scores of right size');
#cmp_ok($hmm->nullTransitions->getdim(1), '==', 2, 'Null transitions of right size');
cmp_ok($hmm->emissions->getdim(0), '==', 260, 'Emissions have right number of states');
cmp_ok($hmm->emissions->getdim(1), '==', 20, 'Emissions have right number of residues');
cmp_ok($hmm->emissions->getdim(2), '==', 2, 'Emissions have both match and insert dimensions');
cmp_ok($hmm->transitions->getdim(0), '==', 260, 'Transitions have right number of states');
cmp_ok($hmm->transitions->getdim(1), '==', 7, 'Transitions have right number of transitions');

cmp_ok($hmm->nullEmissionScores->at(0,0), '==', 2.68618, 'Right null emission score for A');
cmp_ok($hmm->nullEmissionScores->at(0,1), '==', 4.42225, 'Right null emission score for C');
cmp_ok($hmm->nullEmissionScores->at(0,18), '==', 4.58477, 'Right null emission score for W');
cmp_ok($hmm->nullEmissionScores->at(0,19), '==', 3.61503, 'Right null emission score for Y');
is(sprintf('%1.5f', $hmm->nullEmissions->at(0,0)), '0.06814', 'Right null emission prob for A');
is(sprintf('%1.5f', $hmm->nullEmissions->at(0,4)), '0.03132', 'Right null emission prob for F');
is(sprintf('%1.2f', $hmm->nullEmissions->squeeze()->sumover()), '1.00', 'Cross sum of null emissions is 1');

cmp_ok($hmm->emissionScores->at(0, 0, 0), '==', 2.91338, 'Right match emission score for A at 1');
cmp_ok($hmm->emissionScores->at(10, 1, 0), '==', 5.35532, 'Right match emission score for C at 11');
cmp_ok($hmm->emissionScores->at(150, 15, 0), '==', 2.18257, 'Right match emission score for S at 151');
cmp_ok($hmm->emissionScores->at(259, 19, 0), '==', 3.80748, 'Right match emission score for Y at 346');

cmp_ok($hmm->emissionScores->at(0, 0, 1), '==', 2.68618, 'Right insert emission score for A at 1');
cmp_ok($hmm->emissionScores->at(10, 1, 1), '==', 4.42225, 'Right insert emission score for C at 11');
cmp_ok($hmm->emissionScores->at(150, 15, 1), '==', 2.37887, 'Right insert emission score for S at 301');
#cmp_ok($hmm->emissionScores->at(259, 19, 1), '==', -987654321, 'Right insert emission score for Y at 260');

is(sprintf('%1.5f', $hmm->emissions->at(0, 0, 0)), '0.05429', 'Right match emission prob for A at 1');
is(sprintf('%1.5f', $hmm->emissions->at(14, 4, 0)), '0.01433', 'Right match emission prob for F at 15');
is(sprintf('%1.2f', $hmm->emissions->slice('133,:,0')->squeeze()->sumover()), '1.00', 'Cross sum of match emissions at 133 is 1');
is(sprintf('%1.2f', $hmm->emissions->slice('84,:,1')->squeeze()->sumover()), '1.00', 'Cross sum of insert emissions at 84 is 1');