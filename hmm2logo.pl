#!/usr/bin/perl -w
#########
# Author: Benjamin Schuster-Boeckler <bsb@sanger.ac.uk>
# Maintainer: Benjamin Schuster-Boeckler <bsb@sanger.ac.uk>
# Created: ?
# Last Modified: 2004-12-09 (rmp - Sanger port)
#
use strict;

die "Usage: hmm2logo.pl file.hmm pfamid ..." unless @ARGV;
use HMM::Profile;

foreach my $file (@ARGV){
    #($file) = $file =~ /([a-z0-9\/\.\-_]+)/i;

    print "Processing file $file\n";

    my $logo;
    if(-r $file)
    {
        $logo = HMM::Profile->new(-hmmerfile=>$file) || die("Couldn't open $file!\n");
        unless ($file =~ s/h.m$/png/)
	{
	    warn("Could not guess png output name, aborting...\n");
	    exit 1;
	}
    } else {
        print "Starting to parse HMM...\n";
        $logo = HMM::Profile->new(-pfamid=>$file) || die("Couldn't open acc $file!\n");
        print "Finished parsing HMM...\n";
        $file .= '.png';
        print "Name: ".$logo->name."\n";
    }
    my $x_title = 'Relative Entropy';
    my $y_title = 'Contribution';
    my $graph_title = $logo->name();
    my $ysize = 500;
    my $xsize = $logo->length() * 34; #($ysize/12)
    my $greyscale = 0;
    
    print "Drawing Logo...\n";
    $logo->draw_logo(
		     -file	  => $file,
		     -xsize	  => $xsize,
		     -ysize	  => $ysize,
		     -x_title	  => $x_title,
		     -y_title	  => $y_title,
		     -graph_title => $graph_title,
		     -greyscale	  => $greyscale,
		     # This parameter is optional if you have arial installed in
		     # /usr/local/share/fonts/ttfonts
		     -bold_font   => '/Users/schuster/Documents/fonts/arialbd.ttf',
		     -regular_font => '/Users/schuster/Documents/fonts/arial.ttf',
		    ) || die("Error writing $file!\n");
    print "Finished drawing Logo...\n";
}
