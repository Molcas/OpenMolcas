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

subroutine cct3_mc0c1at3b(rowa,cola,rowb,colb,rowc,colc,row,isum,col,a,b,c)
! C = C + A(T)*B

use CCT3_global, only: mhkey
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: rowa, cola, rowb, colb, rowc, colc, row, isum, col
real(kind=wp), intent(in) :: a(rowa,cola), b(rowb,colb)
real(kind=wp), intent(inout) :: c(rowc,colc)
integer(kind=iwp) :: j, k

if (mhkey == 1) then
  ! ESSL
  call DGEMM_('T','N',row,col,isum,One,a,rowa,b,rowb,One,c,rowc)

else
  ! Fortran

  do j=1,col
    do k=1,isum
      c(1:row,j) = c(1:row,j)+a(k,1:row)*b(k,j)
    end do
  end do

end if

return

end subroutine cct3_mc0c1at3b
