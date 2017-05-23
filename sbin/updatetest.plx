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
#***********************************************************************
#
# updatetest.plx (internal use!)
#
# Program that appends generated checkfile to a test input.
# First argument is the test input file and second is the
# checkfile, usage: updatetest.plx inputfile checkfile.
#
# A checkfile is imbedded inside a test file using the EMIL
# file command:
# >>FILE checkfile
# ... checkfile contents ...
# >>EOF
# We just replace the contents with the provided checkfile, or,
# if the ">>FILE checkfile" line is not found it will be added.

die "wrong number of arguments" unless (@ARGV == 2);

my $input = $ARGV[0];
my $check = $ARGV[1];

# read and store the first section of the input

open (INPUT_FH, "<", $input)
  or die "cannot open test input: $!";
my @input_lines;
while (<INPUT_FH>) {
	last if /^>>FILE checkfile$/;
	push @input_lines, $_;
}
close (INPUT_FH);

# write the first section of the input and the contents of the
# checkfile back to the input file.

open (INPUT_FH, ">", $input)
  or die "cannot open test input: $!";
open (CHECK_FH, "<", $check)
  or die "cannot open test input: $!";

foreach my $line (@input_lines) {
	print INPUT_FH $line;
}
print INPUT_FH ">>FILE checkfile\n";
foreach my $line (<CHECK_FH>) {
	print INPUT_FH $line;
}
print INPUT_FH ">>EOF\n";

close (CHECK_FH);
close (INPUT_FH);

exit (0);
