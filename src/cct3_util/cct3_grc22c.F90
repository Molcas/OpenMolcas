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

subroutine cct3_grc22c(a,b,c,mvec,ssa,ssb,pbar,ix)

use CCT3_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: mvec(4096,7), ix
integer(kind=iwp), intent(in) :: ssa, ssb, pbar
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp3, nhelp4, posct, sa1, sa2, sb1, sb2

!1*

if (pbar == 1) then

  ! structure A(p,q)*B(q,r)=C(p,r)

  !1.0 prepare c%d,c%i

  call cct3_grc0(2,0,a%d(0,1),b%d(0,2),0,0,mmul(ssa,ssb),c,posct)

  !1.1 define limitations - no limitations

  !1.2 def symm states and test the limitations

  ix = 0
  do sa1=1,nsym

    sa2 = mmul(ssa,sa1)
    sb1 = sa2

    sb2 = mmul(ssb,sb1)

    !1.3 def mvec,c%d and c%i

    ia = a%i(sa1,1,1)
    ib = b%i(sb1,1,1)
    ic = c%i(sa1,1,1)

    ! yes/no
    if ((a%d(ia,2) > 0) .and. (b%d(ib,2) > 0)) then
      nhelp1 = 1
    else
      cycle
    end if

    ! rowA
    nhelp2 = dimm(a%d(0,1),sa1)

    ! colB
    nhelp3 = dimm(b%d(0,2),sb2)

    ! sum
    nhelp4 = dimm(a%d(0,2),sa2)

    ix = ix+1
    mvec(ix,1) = nhelp1
    mvec(ix,2) = a%d(ia,1)
    mvec(ix,3) = b%d(ib,1)
    mvec(ix,4) = c%d(ic,1)
    mvec(ix,5) = nhelp2
    mvec(ix,6) = nhelp4
    mvec(ix,7) = nhelp3

  end do

end if

return

end subroutine cct3_grc22c
