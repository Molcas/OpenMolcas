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

subroutine Chck_T2p(T21,aSGrp,bSGrp)
! test T2+

implicit none
#include "chcc1.fh"
real*8 T21(1:16*31,1:no*(no+1)/2)
integer aSGrp, bSGrp
integer u, v, a, b, ab, uv, bad, ap, bp
real*8 s

if (aSGrp == 2) then
  ap = nv/2
else
  ap = 0
end if

if (bSGrp == 2) then
  bp = nv/2
else
  bp = 0
end if

bad = 0

uv = 0
do u=1,no
  do v=1,u
    uv = uv+1

    ab = 0
    do a=2,nv/2
      do b=1,a-1
        ab = ab+1

        s = (T2c(bp+b,ap+a,v,u)+T2c(bp+b,ap+a,u,v))/2

        if (abs(T21(ab,uv)-s) > 1.0d-10) bad = bad+1
        T21(ab,uv) = s

      end do
    end do

  end do
end do

if (bad == 0) then
  write(6,*) ' Chck T2+ OK ',bad
else
  write(6,*) ' Chck T2+ Bug !!!!!!! ',bad
end if

return

end subroutine Chck_T2p
