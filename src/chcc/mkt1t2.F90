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

subroutine MkT1T2()
! T1(a,i) = 0 (mozno neskor ine)
! T2(a,b,i,j) = (ai|bj)/Dabij

use chcc_global, only: no, nv, OEo, OEv, Q21, T1c, T2c
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: a, b, i, j

call mma_allocate(T1c,nv,no,label='T1c')
call mma_allocate(T2c,nv,nv,no,no,label='T2c')

do i=1,no
  do a=1,nv
    T1c(a,i) = Zero
  end do
end do

do j=1,no
  do i=1,no
    do b=1,nv
      do a=1,nv
        T2c(a,b,i,j) = Q21(a,i,b,j)/(OEo(i)+OEo(j)-OEv(a)-OEv(b))
      end do
    end do
  end do
end do

return

end subroutine MkT1T2
