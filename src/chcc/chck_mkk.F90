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

subroutine Chck_mkK()
! make K(be,u,i,a)

use chcc_global, only: Kc, no, nv, Q1, Q21, Q22, Q3, T1c, T2c
use stdalloc, only: mma_allocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: a, b, be, i, j, u
real(kind=wp) :: s

call mma_allocate(Kc,no,nv,no,nv,label='Kc')

do a=1,nv
  do i=1,no
    do u=1,no
      do be=1,nv

        s = Zero

        s = s+Q22(be,a,i,u)

        do j=1,no
          s = s-Q1(a,j,i,u)*T1c(be,j)
        end do

        do b=1,nv
          s = s+Q3(a,be,b,i)*T1c(b,u)
        end do

        do j=1,no
          do b=1,nv
            s = s-Q21(b,i,a,j)*(T2c(b,be,u,j)*Half+T1c(b,u)*T1c(be,j))
          end do
        end do

        Kc(i,be,u,a) = s

      end do
    end do
  end do
end do

write(u6,*) ' K done '

return

end subroutine Chck_mkK
