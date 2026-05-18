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

subroutine PCOLLVEC(IVEC,iTYPE)

use definitions, only: iwp
use caspt2_module, only: nCases, nSym, nInDep, nASup, nISup

implicit none
integer(kind=iwp), intent(in) :: iVec, iType
integer(kind=iwp) iCase, iSym, NAS, NIS, NW

!***********************************************************************
do ICASE=1,NCASES
  do ISYM=1,NSYM
    if (NINDEP(ISYM,ICASE) == 0) cycle
    if (ITYPE == 0) then
      NAS = NINDEP(ISYM,ICASE)
    else
      NAS = NASUP(ISYM,ICASE)
    end if
    NIS = NISUP(ISYM,ICASE)
    NW = NAS*NIS
    if (NW == 0) cycle
    call DRA2SOLV(NAS,NIS,iCASE,iSYM,iVEC)
  end do
end do

end subroutine PCOLLVEC
