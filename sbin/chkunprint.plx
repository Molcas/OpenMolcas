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
# filter to check unprintable chars in a file
#
# returns 1 if there are some Ctrl chars, or several extended ASCII

$thename=$ARGV[0] || '';
help () if ($thename eq '-h');
$shownum=0;
if($thename eq '-n') {$shownum=1;}
$i=0;
while(<STDIN>)
{
  $i++;
  if($_=~/Segmentation *Vio/)
  {
    print "Segmentation Error\n";
    exit(2);
  }

  if ($_=~/[\x00-\x08\x0e-\x1F]/ )
  {
    print "garbage in line= $i\n" if($shownum==1);
    exit(1);
  }

  if ( $_=~/[\xA0-\xFF].*[\xA0-\xFF].*[\xA0-\xFF]/)
  {
    print "garbage in line= $i\n" if($shownum==1);
    exit(1);
  }

}

exit (0);
sub help {
  print "Usage: molcas chkunprint.plx <file\n\n";

  print "returns 0 if there is no garbage in the file,\n";
  print "     or 1\n";
  print " To print line number use flag -n\n";
  exit(0);
}
