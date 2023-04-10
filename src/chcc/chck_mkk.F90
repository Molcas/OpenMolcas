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

implicit none
#include "chcc1.fh"
integer be, u, i, a
integer j, b
real*8 s

do a=1,nv
  do i=1,no
    do u=1,no
      do be=1,nv

        s = 0.0d0

        s = s+Q22(be,a,i,u)

        do j=1,no
          s = s-Q1(a,j,i,u)*T1c(be,j)
        end do

        do b=1,nv
          s = s+Q3(a,be,b,i)*T1c(b,u)
        end do

        do j=1,no
          do b=1,nv
            s = s-Q21(b,i,a,j)*(T2c(b,be,u,j)/2+T1c(b,u)*T1c(be,j))
          end do
        end do

        Kc(i,be,u,a) = s

      end do
    end do
  end do
end do

write(6,*) ' K done '

return

end subroutine Chck_mkK
