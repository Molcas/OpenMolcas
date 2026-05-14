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

subroutine PRSM2(A,NDIM)
! PRINT LOWER TRIANGULAR MATRIX PACKED IN COLUMN WISE FASHION

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NDIM
real(kind=wp), intent(in) :: A(nTri_Elem(NDIM))
integer(kind=iwp) :: I, J

do I=1,NDIM
  write(u6,1010) I,(A((J-1)*NDIM-nTri_Elem(J-1)+I),J=1,I)
end do

return
1010 format('0',2X,I3,5(ES14.7),/,(1X,5X,5(ES14.7)))

end subroutine PRSM2
