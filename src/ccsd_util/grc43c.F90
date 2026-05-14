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

subroutine grc43C(a,b,c,mvec,ssa,ssb,pbar,ix)

use ccsd_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: mvec(4096,7)
integer(kind=iwp), intent(in) :: ssa, ssb, pbar
integer(kind=iwp), intent(inout) :: ix
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp21, nhelp22, nhelp3, nhelp4, nhelp41, nhelp42, nsyma2, ntest1, ntest2, &
                     posct, sa1, sa12, sa123, sa2, sa3, sa4, sb1, sb12, sb2, sb3

!1*

if (pbar == 1) then

  ! structure A(p,qrs)*B(qrs)=C(p)
  ! implemented in grc43y

else if (pbar == 2) then

  ! structure A(pq,rs)*B(rs,t)=C(pq,t)

  !1.1 define limitations - p>q,r,s must be tested - ntest1
  !                       - p,q,r>s must be tested - ntest2

  if ((a%d(0,6) == 1) .or. (a%d(0,6) == 4)) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if ((a%d(0,6) == 3) .or. (a%d(0,6) == 4)) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

  !1.0 prepare c%d,c%i

  call grc0(3,ntest1,a%d(0,1),a%d(0,2),b%d(0,3),0,mmul(ssa,ssb),posct,c)

  !1.2 def symm states and test the limitations

  ix = 1
  do sa1=1,nsym
    if (ntest1 == 1) then
      nsyma2 = sa1
    else
      nsyma2 = nsym
    end if

    do sa2=1,nsyma2
      sa12 = mmul(sa1,sa2)

      do sa3=1,nsym
        sa123 = mmul(sa12,sa3)
        sb1 = sa3

        sa4 = mmul(ssa,sa123)
        sb2 = sa4
        sb12 = mmul(sb1,sb2)
        ! Meggie out
        if ((ntest2 == 1) .and. (sa3 < sa4)) cycle

        sb3 = mmul(ssb,sb12)

        !1.3 def mvec,c%d and c%i

        ia = a%i(sa1,sa2,sa3)
        ib = b%i(sb1,sb2,1)
        ic = c%i(sa1,sa2,1)

        ! yes/no
        if ((a%d(ia,2) <= 0) .or. (b%d(ib,2) <= 0)) cycle
        nhelp1 = 1

        ! rowA
        nhelp21 = dimm(a%d(0,1),sa1)
        nhelp22 = dimm(a%d(0,2),sa2)
        if ((ntest1 == 1) .and. (sa1 == sa2)) then
          nhelp2 = nhelp21*(nhelp21-1)/2
        else
          nhelp2 = nhelp21*nhelp22
        end if

        ! colB
        nhelp3 = dimm(b%d(0,3),sb3)

        ! sum
        nhelp41 = dimm(a%d(0,2),sa2)
        nhelp42 = dimm(a%d(0,3),sa3)
        if ((ntest2 == 1) .and. (sa2 == sa3)) then
          nhelp4 = nhelp41*(nhelp41-1)
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

end if
ix = ix-1

return

end subroutine grc43C
