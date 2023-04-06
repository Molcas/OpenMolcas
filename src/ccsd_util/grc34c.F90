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

subroutine grc34C(a,b,c,mvec,ssa,ssb,pbar,ix)

use ccsd_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: mvec(4096,7)
integer(kind=iwp), intent(in) :: ssa, ssb, pbar
integer(kind=iwp), intent(inout) :: ix
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp3, nhelp31, nhelp32, nhelp4, nhelp41, nhelp42, ntest1, ntest2, posct, sa1, &
                     sa12, sa2, sa3, sb1, sb12, sb123, sb2, sb3, sb4

!1*

if (pbar == 1) then

  ! structure A(p,qr)*B(qr,st)=C(p,st)

  !1.1 define limitations - q>r,s,t must be tested - ntest1
  !                       - q,r,s>t must be tested - ntest2

  if ((b%d(0,6) == 1) .or. (b%d(0,6) == 4)) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if ((b%d(0,6) == 3) .or. (b%d(0,6) == 4)) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

  !1.0 prepare c%d,c%i

  if (ntest2 == 1) then
    nhelp1 = 2
  else
    nhelp1 = 0
  end if

  call grc0(3,nhelp1,a%d(0,1),b%d(0,3),b%d(0,4),0,mmul(ssa,ssb),posct,c)

  !1.2 def symm states and test the limitations

  ix = 1
  do sa1=1,nsym

    do sa2=1,nsym
      sa12 = mmul(sa1,sa2)
      sb1 = sa2

      sa3 = mmul(ssa,sa12)
      sb2 = sa3
      sb12 = mmul(sb1,sb2)
      ! Meggie out
      if ((ntest1 == 1) .and. (sa2 < sa3)) cycle

      do sb3=1,nsym
        sb123 = mmul(sb12,sb3)

        sb4 = mmul(ssb,sb123)
        ! Meggie out
        if ((ntest2 == 1) .and. (sb3 < sb4)) cycle

        !1.3 def mvec,c%d and c%i

        ia = a%i(sa1,sa2,1)
        ib = b%i(sb1,sb2,sb3)
        ic = c%i(sa1,sb3,1)

        ! yes/no
        if ((a%d(ia,2) <= 0) .or. (b%d(ib,2) <= 0)) cycle
        nhelp1 = 1

        ! rowA
        nhelp2 = dimm(a%d(0,1),sa1)

        ! colB
        nhelp31 = dimm(b%d(0,3),sb3)
        nhelp32 = dimm(b%d(0,4),sb4)
        if ((ntest2 == 1) .and. (sb3 == sb4)) then
          nhelp3 = nhelp31*(nhelp31-1)/2
        else
          nhelp3 = nhelp31*nhelp32
        end if

        ! sum
        nhelp41 = dimm(a%d(0,2),sa2)
        nhelp42 = dimm(a%d(0,3),sa3)
        if ((ntest1 == 1) .and. (sa2 == sa3)) then
          nhelp4 = nhelp41*(nhelp41-1)/2
        else
          nhelp4 = nhelp41*nhelp42
        end if

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

else if (pbar == 3) then

! structure A(pq,r)*B(r,stu)=C(pq,stu)
! not used

end if
ix = ix-1

return

end subroutine grc34C
