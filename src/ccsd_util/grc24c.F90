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

subroutine grc24C(a,b,c,mvec,ssa,ssb,pbar,ix)

use ccsd_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: mvec(4096,7)
integer(kind=iwp), intent(in) :: ssa, ssb, pbar
integer(kind=iwp), intent(inout) :: ix
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp3, nhelp31, nhelp32, nhelp33, nhelp4, nsymb3, ntest1, ntest2, posct, sa1, &
                     sa2, sb1, sb12, sb123, sb2, sb3, sb4

!1*

if (pbar == 1) then

  ! structure A(p,q)*B(q,rst)=C(p,rst)

  !1.0 prepare c%d,c%i

  call grc0(4,b%d(0,6),a%d(0,1),b%d(0,2),b%d(0,3),b%d(0,4),mmul(ssa,ssb),posct,c)

  !1.1 define limitations - q,r>s,t must be tested - ntest1
  !                       - q,r,s>t must be tested - ntest2

  if (b%d(0,6) == 2) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if (b%d(0,6) == 3) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

  !1.2 def symm states and test the limitations

  ix = 1
  do sa1=1,nsym

    sa2 = mmul(ssa,sa1)
    sb1 = sa2

    do sb2=1,nsym
      sb12 = mmul(sb1,sb2)
      if (ntest1 == 1) then
        nsymb3 = sb2
      else
        nsymb3 = nsym
      end if

      do sb3=1,nsymb3
        sb123 = mmul(sb12,sb3)

        sb4 = mmul(ssb,sb123)
        ! Meggie out
        if ((ntest2 == 1) .and. (sb3 < sb4)) cycle

        !1.3 def mvec,c%d and c%i

        ia = a%i(sa1,1,1)
        ib = b%i(sb1,sb2,sb3)
        ic = c%i(sa1,sb2,sb3)

        ! yes/no
        if ((a%d(ia,2) <= 0) .or. (b%d(ib,2) <= 0)) cycle
        nhelp1 = 1

        ! rowA
        nhelp2 = dimm(a%d(0,1),sa1)

        ! colB
        nhelp31 = dimm(b%d(0,2),sb2)
        nhelp32 = dimm(b%d(0,3),sb3)
        nhelp33 = dimm(b%d(0,4),sb4)
        if ((ntest1 == 1) .and. (sb2 == sb3)) then
          nhelp3 = nhelp31*(nhelp31-1)*nhelp33/2
        else if ((ntest2 == 1) .and. (sb3 == sb4)) then
          nhelp3 = nhelp31*nhelp32*(nhelp32-1)/2
        else
          nhelp3 = nhelp31*nhelp32*nhelp33
        end if

        ! sum
        nhelp4 = dimm(a%d(0,2),sa2)

        mvec(ix,1) = nhelp1
        mvec(ix,2) = a%d(ia,1)
        mvec(ix,3) = b%d(ib,1)
        mvec(ix,4) = c%d(ic,1)
        mvec(ix,5) = nhelp2
        mvec(ix,6) = nhelp4
        mvec(ix,7) = nhelp3

        ix = ix+1

      end do
    end do
  end do

end if
ix = ix-1

return

end subroutine grc24C
