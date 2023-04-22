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

subroutine Chck_T21od(T21,beSGrp,gaSGrp)
! test T2n+

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, nv, Q4, T2c
use Constants, only: Zero, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: beSGrp, gaSGrp
real(kind=wp), intent(inout) :: T21(nv/2,nv/2,nTri_Elem(no-1))
integer(kind=iwp) :: a, b, bad, be, bega, bep, ga, gap, u, uv, v
real(kind=wp) :: s

if (beSGrp == 2) then
  bep = nv/2
else
  bep = 0
end if

if (gaSGrp == 2) then
  gap = nv/2
else
  gap = 0
end if

bad = 0

uv = 0
do u=2,no
  do v=1,u-1
    uv = uv+1

    bega = 0
    do be=1,nv/2
      do ga=1,nv/2
        bega = bega+1

        s = Zero
        do a=1,nv
          b = a
          s = s+(Q4(b,gap+ga,a,bep+be)+Q4(b,bep+be,a,gap+ga))*(T2c(b,a,v,u)+T2c(b,a,u,v))*Quart
        end do

        s = Zero
        do a=2,nv
          do b=1,a-1
            s = s+(Q4(b,gap+ga,a,bep+be)-Q4(b,bep+be,a,gap+ga))*(T2c(b,a,v,u)-T2c(b,a,u,v))*Half
          end do
        end do

        if (abs(T21(be,ga,uv)-s) > 1.0e-10_wp) then
          bad = bad+1
          !write(u6,99) be,ga,u,v
          !99 format(4(i3,1x))
        end if
        T21(be,ga,uv) = s

      end do
    end do

  end do
end do

if (bad == 0) then
  write(u6,*) ' Chck T2 OK ',bad
else
  write(u6,*) ' Chck T2 Bug !!!!!!! ',bad
end if

return

end subroutine Chck_T21od
