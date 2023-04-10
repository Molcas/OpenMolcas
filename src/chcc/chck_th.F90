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

subroutine Chck_Th(T2)
! check T2 = T2(a,b,j_b,u_a)

implicit none
#include "chcc1.fh"
real*8 T2(1:nv*(nv+1)/2,1:no,1:no)
! help var
integer u, a, b, j, bad, ab
real*8 s

bad = 0

do u=1,no
  do j=1,no
    ab = 0
    do a=1,nv
      do b=1,a
        ab = ab+1

        s = T2c(a,b,j,u)+T1c(a,j)*T1c(b,u)

        if (abs(T2(ab,j,u)-s) > 1.0d-10) then
          T2(ab,j,u) = s
          bad = bad+1
        end if

      end do
    end do
  end do
end do

write(6,*) ' T2  Chck :',bad

return

end subroutine Chck_Th
