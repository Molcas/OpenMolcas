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

subroutine mc0c1at3b(rowa,cola,rowb,colb,rowc,colc,row,rsum,col,a,b,c)
! C= A(T)*B

#include "ccsd1.fh"
integer rowa, cola, rowb, colb, rowc, colc
integer row, rsum, col
real*8 a(1:rowa,1:cola)
real*8 b(1:rowb,1:colb)
real*8 c(1:rowc,1:colc)
! help variables
integer i, j, k

if (mhkey == 1) then
  ! ESSL
  call DGEMM_('T','N',row,col,rsum,1.0d0,a,rowa,b,rowb,1.0d0,c,rowc)

else
  ! Fortran

  do j=1,col
    do i=1,row
      do k=1,rsum
        c(i,j) = c(i,j)+a(k,i)*b(k,j)
      end do
    end do
  end do

end if

return

end subroutine mc0c1at3b
