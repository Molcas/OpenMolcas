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

subroutine Chck_W1(W1,aSGrp,beSGrp,bSGrp,gaSGrp)
! cek W1

use Index_Functions, only: nTri_Elem
use chcc_global, only: nv, Q4
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: W1(nTri_Elem(nv/2-1),nTri_Elem(nv/2))
integer(kind=iwp), intent(in) :: aSGrp, beSGrp, bSGrp, gaSGrp
integer(kind=iwp) :: a, ab, ap, b, bad, be, bega, bep, bp, ga, gap
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

if (gaSGrp == 2) then
  gap = nv/2
else
  gap = 0
end if

if (beSGrp == 2) then
  bep = nv/2
else
  bep = 0
end if

bad = 0
bega = 0
do be=1,nv/2
  do ga=1,be
    bega = bega+1
    ab = 0
    do a=2,nv/2
      do b=1,a-1
        ab = ab+1
        s = (Q4(ap+a,bep+be,bp+b,gap+ga)+Q4(ap+a,gap+ga,bp+b,bep+be))/1
        if (abs(W1(ab,bega)-s) > 1.0e-10_wp) then
          bad = bad+1
          !write(u6,99) a,b,be,ga,ab,bega,s,W1(a,be,ga)
          !99 format(4(i2,1x),2(i6,1x),2(f15.10))
        end if
        W1(ab,bega) = s
      end do
    end do
  end do
end do

if (bad == 0) then
  write(u6,*) ' Chck W OK ',bad
else
  write(u6,*) ' Chck W Bug !!!!!!! ',bad
end if

return

end subroutine Chck_W1
