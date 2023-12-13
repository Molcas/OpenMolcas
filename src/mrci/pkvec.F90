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

subroutine PKVEC(NITEM,CVEC,ICVEC)
!***********************************************************************
!
! PURPOSE:
! ENCODE THE CI-VECTOR BY CHANGING THE NUMBER REPRESENTATION FROM
! REAL*8 ==> INTEGER*4
!
! NOTE:
! THE INCOMING DATA CVEC SHOULD NOT BE GREATER THAN 1.0.
! THE ACCURACY OF THE UNPACKED VALUES IS APPROX. 1.0e-9.
!
!**** M.P. FUELSCHER, UNIVERSITY OF LUND, SWEDEN, NOV. 1990 ************

use Definitions, only: wp, iwp
#ifdef _CYGWIN_
use Definitions, only: r4
#endif

implicit none
integer(kind=iwp), intent(in) :: NITEM
real(kind=wp), intent(in) :: CVEC(NITEM)
integer(kind=iwp), intent(out) :: ICVEC(NITEM)
integer(kind=iwp) :: ITEM
real(kind=wp), parameter :: SCL = 2147483647.0_wp
#ifdef _CYGWIN_
! NOTE VERY CAREFULLY!! NINT() (and maybe similar) intrinsics
! are BROKEN in CYGWIN GFORTAN!!
! So the intermediate variable tmp is necessary for it to work!
real(kind=r4) :: tmp
#endif

do ITEM=1,NITEM
# ifdef _CYGWIN_
  tmp = SCL*CVEC(ITEM)
  ICVEC(ITEM) = nint(tmp)
# else
  ICVEC(ITEM) = nint(SCL*CVEC(ITEM))
# endif
end do

return

end subroutine PKVEC
