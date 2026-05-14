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

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, nv, T2c
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: T21(nTri_Elem(nv/2-1),nTri_Elem(no))
integer(kind=iwp), intent(in) :: aSGrp, bSGrp
integer(kind=iwp) :: a, ab, ap, b, bad, bp, u, uv, v
real(kind=wp) :: s

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

        s = (T2c(bp+b,ap+a,v,u)+T2c(bp+b,ap+a,u,v))*Half

        if (abs(T21(ab,uv)-s) > 1.0e-10_wp) bad = bad+1
        T21(ab,uv) = s

      end do
    end do

  end do
end do

if (bad == 0) then
  write(u6,*) ' Chck T2+ OK ',bad
else
  write(u6,*) ' Chck T2+ Bug !!!!!!! ',bad
end if

return

end subroutine Chck_T2p
