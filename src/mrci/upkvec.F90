!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990, Markus P. Fuelscher                              *
!***********************************************************************

subroutine UPKVEC(NITEM,ICVEC,CVEC)
!***********************************************************************
!
! PURPOSE:
! DECODE THE CI-VECTOR BY CHANGING THE NUMBER REPRESENTATION FROM
! INTEGER*4 ==> REAL*8
!
!**** M.P. FUELSCHER, UNIVERSITY OF LUND, SWEDEN, NOV. 1990 ************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NITEM, ICVEC(NITEM)
real(kind=wp), intent(out) :: CVEC(NITEM)
integer(kind=iwp) :: ITEM
real(kind=wp), parameter :: SCL = One/2147483647.0_wp

do ITEM=1,NITEM
  CVEC(ITEM) = SCL*real(ICVEC(ITEM),kind=wp)
end do

return

end subroutine UPKVEC
