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

subroutine Chck_Tx(T)
! check T(a,b,i,j)

implicit none
#include "chcc1.fh"
!real*8 T(1:nv,1:nv,1:no,1:no)
real*8 T(1:nv,1:no,1:nv,1:no)
!real*8 T(1:nv,1:nv,1:no)
integer b, j, a, i, bad
real*8 s

bad = 0
do j=1,no
  do i=1,no
    do b=1,nv
      do a=1,nv

        s = T2c(a,b,i,j)

        if (abs(T(b,i,a,j)-s) > 1.0d-10) bad = bad+1

      end do
    end do
  end do
end do

write(6,*) ' Chck T2 :',bad

return

end subroutine Chck_Tx
