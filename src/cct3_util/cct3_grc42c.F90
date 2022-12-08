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

subroutine cct3_grc42c(a,b,c,mvec,ssa,ssb,pbar,ix)

use CCT3_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: mvec(4096,7), ix
integer(kind=iwp), intent(in) :: ssa, ssb, pbar
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp21, nhelp22, nhelp23, nhelp3, nhelp4, nsyma2, nsyma3, ntest1, ntest2, posct, &
                     sa1, sa12, sa123, sa2, sa3, sa4, sb1, sb2

!1*

if (pbar == 1) then

  ! not possible

else if (pbar == 2) then

  ! structure A(pq,rs)*B(rs)=C(pq)
  ! implemented in grc42y

else if (pbar == 3) then

  ! structure A(pqr,s)*B(s,t)=C(pqr,t)

  !1.1 define limitations -  p>q,r,s must be tested - ntest1
  !                        -  p,q>r,s must be tested - ntest2

  !1.0 prepare c%d,c%i

  call cct3_grc0(4,a%d(0,6),a%d(0,1),a%d(0,2),a%d(0,3),b%d(0,2),mmul(ssa,ssb),c,posct)

  if (a%d(0,6) == 1) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if (a%d(0,6) == 2) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

  !1.2 def symm states and test the limitations

  ix = 0
  do sa1=1,nsym
    if (ntest1 == 1) then
      nsyma2 = sa1
    else
      nsyma2 = nsym
    end if

    do sa2=1,nsyma2
      sa12 = mmul(sa1,sa2)
      if (ntest2 == 1) then
        nsyma3 = sa2
      else
        nsyma3 = nsym
      end if

      do sa3=1,nsyma3
        sa123 = mmul(sa12,sa3)

        sa4 = mmul(ssa,sa123)
        sb1 = sa4

        sb2 = mmul(ssb,sb1)

        !1.3 def mvec,c%d and c%i

        ia = a%i(sa1,sa2,sa3)
        ib = b%i(sb1,1,1)
        ic = c%i(sa1,sa2,sa3)

        ! yes/no
        if ((a%d(ia,2) > 0) .and. (b%d(ib,2) > 0)) then
          nhelp1 = 1
        else
          cycle
        end if

        ! rowA
        nhelp21 = dimm(a%d(0,1),sa1)
        nhelp22 = dimm(a%d(0,2),sa2)
        nhelp23 = dimm(a%d(0,3),sa3)
        if ((ntest1 == 1) .and. (sa1 == sa2)) then
          nhelp2 = nhelp21*(nhelp21-1)*nhelp23/2
        else if ((ntest2 == 1) .and. (sa2 == sa3)) then
          nhelp2 = nhelp21*nhelp22*(nhelp22-1)/2
        else
          nhelp2 = nhelp21*nhelp22*nhelp23
        end if

        ! colB
        nhelp3 = dimm(b%d(0,2),sb2)

        ! sum
        nhelp4 = dimm(a%d(0,4),sa4)

        ix = ix+1
        mvec(ix,1) = nhelp1
        mvec(ix,2) = a%d(ia,1)
        mvec(ix,3) = b%d(ib,1)
        mvec(ix,4) = c%d(ic,1)
        mvec(ix,5) = nhelp2
        mvec(ix,6) = nhelp4
        mvec(ix,7) = nhelp3

      end do
    end do
  end do

end if

return

end subroutine cct3_grc42c
