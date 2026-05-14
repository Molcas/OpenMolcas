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

subroutine Chck_A(AA)
! check AA(ij,u,v)

use Index_Functions, only: nTri_Elem
use chcc_global, only: Ac, no, nv, Q0, Q1, Q21, T1c, T2c
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: AA(nTri_Elem(no),no,no)
integer(kind=iwp) :: a, b, bad, i, ij, j, u, v
real(kind=wp) :: s

call mma_allocate(Ac,no,no,no,no,label='Ac')

bad = 0

ij = 0
do i=1,no
  do j=1,i
    ij = ij+1

    do u=1,no
      do v=1,no

        s = Q0(u,i,v,j)

        do a=1,nv
          s = s+Q1(a,j,u,i)*T1c(a,v)
        end do

        do a=1,nv
          s = s+Q1(a,i,v,j)*T1c(a,u)
        end do

        do a=1,nv
          do b=1,nv
            s = s+Q21(a,i,b,j)*(T2c(a,b,u,v)+T1c(a,u)*T1c(b,v))
          end do
        end do

        Ac(i,j,u,v) = s

        !write(u6,99) i,j,u,v,AA(ij,u,v),s,AA(ij,u,v)-s
        !99 format(4(i2,1x),3(f15.10,1x))
        if (abs(AA(ij,u,v)-s) > 1.0e-10_wp) bad = bad+1

      end do
    end do
  end do
end do

do i=2,no
  Ac(1:i-1,i,:,:) = Ac(i,1:i-1,:,:)
end do

write(u6,*) ' A   Chck :',bad

return

end subroutine Chck_A
