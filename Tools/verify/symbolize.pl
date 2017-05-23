#!/usr/bin/perl
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

# script to convert stack trace with pc addresses
# to source code files and line numbers

print "symbolized stack trace:\n";
while (<STDIN>) {
    if (/^(\s+#\d+)\s+0x[0-9a-f]+\s+\(([^\(\)]*)\)/) {
        my $counter = $1;
        my $address = $2;
        my ($exe,$pc) = split(/\+/, $address);
        my $source = `addr2line -f -p -e $exe $pc`;
        if ($source =~ /\?\?:\?/) {
            print "$counter (?) $exe $pc\n";
        } else {
            print "$counter $source";
        }
    }
}
