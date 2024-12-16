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

module splitcas_data

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iDimBlockA, iDimBlockACNF, iterSplit, lrootSplit, MxIterSplit
real(kind=wp) :: EnInSplit, GapSpli, PercSpli, thrSplit
logical(kind=iwp) :: DoSplitCAS, EnerSplit, FOrdSplit, NumSplit, PerSplit

public :: DoSplitCAS, EnerSplit, EnInSplit, FordSplit, gapSpli, iDimBlockA, iDimBlockACNF, iterSplit, lRootSplit, MxIterSplit, &
          NumSplit, percSpli, PerSplit, ThrSplit

end module splitcas_data
