require 5.006;

use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'           => 'HMM',
    'VERSION_FROM'   => 'HMM.pm', # finds $VERSION
    'PERL_MALLOC_OK' => 1, # should be able to use perl malloc()
    'PREREQ_PM'      => { PDL             => '2.3.4',
                          Imager          => '0.4.3',
                          LWP             => '5.10',
                          HTTP::Request   => '1.0',
                          HTTP::Response  => '1.0',
                          Test::Simple    => '0.86' },
    'dist'      => {COMPRESS     => 'gzip -9f',
                    SUFFIX       => '.gz', 
                    DIST_DEFAULT => 'all tardist'}
);
