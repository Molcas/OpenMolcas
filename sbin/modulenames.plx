#!/usr/bin/env perl
#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2014, Steven Vancoillie                                *
#***********************************************************************
#
# modulenames.plx:
#
# extract module names from a test input and generate a proper test name for
# use with CMake (no spaces, special chars)
#
# Steven Vancoillie, Lund, 20th of June 2014

# Perl modules
use warnings;
use strict;

sub usage {
    print <<USAGE;

 Usage: modulenames.plx <filename>

        Program to generate a test name by concatenating the
        module names used in a test input given by filename.

USAGE
    exit 1;
}

# only argument is the file name
my $filename = $ARGV[0];
&usage unless ($filename);

# open the file
my %module_names;
open TESTINPUT, '<', $filename or die "Error: failed to open file $filename\n";
while (<TESTINPUT>) {
        my @matches = ($_ =~ /&(\w+)/g);
        foreach my $name (@matches) {
                my $lc_name = lc $name;
                next if $lc_name =~ /end/;
                $module_names{$lc_name}++;
        }
}
close TESTINPUT;

my $output = join ('_', keys %module_names);

print "$output\n";

exit 0;
