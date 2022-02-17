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

subroutine CSCALE(INDX,INTSYM,C,X)

use mrci_global, only: IRC, LSYM, NDIAG, NVIRT
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: INDX(*), INTSYM(*)
real(kind=wp), intent(inout) :: C(*)
real(kind=wp), intent(in) :: X
integer(kind=iwp) :: II1, MA, NA
integer(kind=iwp), external :: JSUNP

do II1=IRC(3)+1,IRC(4)
  if (JSUNP(INTSYM,II1) == LSYM) then
    NA = INDX(II1)
    do MA=1,NVIRT
      C(NA+NDIAG(MA)) = X*C(NA+NDIAG(MA))
    end do
  end if
end do

return

end subroutine CSCALE
