#!/usr/bin/env perl

use strict;
use warnings;
use HMM::Profile;
use Data::Printer;

my $file = shift;
my $height_logodds = shift || 0;

my $logo = HMM::Profile->new( -hmmerfile => $file ) or 
   die "Failed in making HMM logo from $file!\n";

my $outfile     = "$file.png"; #"hmmLogo.png";
my $graph_title = $logo->name();
my $ysize       = 500;
my $xsize       = $logo->length() * 34;
my $greyscale   = 0;

if (0==1) {
$logo->print_logo_dimensions(
    -xsize          => $xsize,
    -ysize          => $ysize,
    -x_title        => 'Relative Entropy',
    -y_title        => 'Contribution',
    -graph_title    => $graph_title,
    -greyscale      => $greyscale,
    -height_logodds => $height_logodds
  )  or die "Error writing $outfile!\n";
}

#Now go and make the logos
print STDOUT "Drawing Logo...\n";
$logo->draw_logo(
    -file           => $outfile,
    -xsize          => $xsize,
    -ysize          => $ysize,
    -x_title        => 'Relative Entropy',
    -y_title        => 'Contribution',
    -graph_title    => $graph_title,
    -greyscale      => $greyscale,
    -height_logodds => $height_logodds
  )  or die "Error writing $outfile!\n";

my $data = $logo->flat($height_logodds);
#print STDOUT p( $data);
print STDOUT "Finished drawing Logo...\n";
