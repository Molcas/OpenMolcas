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

integer function Cho_P_LocalSP(iShlAB)
!
! Purpose: return local shell pair corresponding to global shell
!          pair iShlAB (returns 0 if not found).

use ChoArr, only: MySP, n_MySP

implicit none
integer iShlAB
#include "cho_para_info.fh"
integer iSP

if (Cho_Real_Par) then
  Cho_P_LocalSP = 0
  do iSP=1,n_mySP
    if (mySP(iSP) == iShlAB) then
      Cho_P_LocalSP = iSP
      return
    end if
  end do
else
  Cho_P_LocalSP = iShlAB
end if

end function Cho_P_LocalSP
