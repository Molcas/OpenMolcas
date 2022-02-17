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

subroutine COUNT_MRCI(NINTGR,NSYM,NORB,MUL)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: NINTGR
integer(kind=iwp), intent(in) :: NSYM, NORB(*), MUL(8,8)
integer(kind=iwp) :: IJS, IS, ISUM, JS, NDPROD(8), NORBT

! COUNT TWO-ELECTRON INTEGRALS
! FIRST, COUNT NUMBER OF DENSITY PRODUCTS IN EACH SYMMETRY:
NDPROD(1) = 0
NORBT = 0
do IS=1,NSYM
  NDPROD(IS) = 0
  NORBT = NORBT+NORB(IS)
end do
do IJS=1,NSYM
  ISUM = 0
  do IS=1,NSYM
    JS = MUL(IS,IJS)
    if (JS <= IS) ISUM = ISUM+NORB(IS)*NORB(JS)
  end do
  NDPROD(IJS) = ISUM
end do
NDPROD(1) = (NDPROD(1)+NORBT)/2
! THEN COUNT NUMBER OF TOTALLY SYMMETRIC PRODUCTS OF DENS-PRODUCTS:
NINTGR = 0
do IJS=1,NSYM
  NINTGR = NINTGR+NDPROD(IJS)*(1+NDPROD(IJS))
end do
NINTGR = NINTGR/2

return

end subroutine COUNT_MRCI
