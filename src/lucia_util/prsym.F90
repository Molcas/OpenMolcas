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

subroutine PRSYM(A,MATDIM)
! PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
! THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MATDIM
real(kind=wp), intent(in) :: A(nTri_Elem(MATDIM))
integer(kind=iwp) :: I, J, JSTART, JSTOP

JSTART = 1
JSTOP = 0
do I=1,MATDIM
  JSTART = JSTART+I-1
  JSTOP = JSTOP+I
  write(u6,1010) I,(A(J),J=JSTART,JSTOP)
end do

return
1010 format('0',2X,I3,5(ES14.7),/,(1X,5X,5(ES14.7)))

end subroutine PRSYM
