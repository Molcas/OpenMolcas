!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! Stuff from splitcas.fh
module splitcas_data
implicit none
Logical DoSplitCAS, EnerSplit, NumSplit, PerSplit, FOrdSplit
integer lrootSplit, MxIterSplit
integer iDimBlockA, iDimBlockACNF, iterSplit
real*8 thrSplit, GapSpli, PercSpli
real*8 EnInSplit
Integer, Parameter :: mxDimBlockA = 2000
Real*8, Parameter :: min_ThrSplit = 1.0d-12
end module splitcas_data
