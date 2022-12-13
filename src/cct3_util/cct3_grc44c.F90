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

subroutine cct3_grc44c(a,b,c,mvec,ssa,ssb,pbar,ix)

use CCT3_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: mvec(4096,7), ix
integer(kind=iwp), intent(in) :: ssa, ssb, pbar
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp21, nhelp22, nhelp3, nhelp31, nhelp32, nhelp4, nhelp41, nhelp42, nhelp43, &
                     nsyma2, nsyma3, ntest1, ntest2, ntest3, posct, sa1, sa12, sa123, sa2, sa3, sa4, sb1, sb12, sb123, sb2, sb3, sb4

!1*

if (pbar == 1) then

  ! structure A(p,qrs)*B(qrs,t)=C(p,t)

  !1.0 prepare c%d,c%i

  call cct3_grc0(2,0,a%d(0,1),b%d(0,4),0,0,mmul(ssa,ssb),c,posct)

  !1.1 define limitations - p,q>r,s - ntest1, p,q,r>s - ntest2

  if (a%d(0,6) == 2) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if (a%d(0,6) == 3) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

  !1.2 def symm states and test the limitations

  ix = 0
  do sa1=1,nsym

    do sa2=1,nsym
      sa12 = mmul(sa1,sa2)
      sb1 = sa2
      if (ntest1 == 1) then
        nsyma3 = sa2
      else
        nsyma3 = nsym
      end if

      do sa3=1,nsyma3
        sa123 = mmul(sa12,sa3)
        sb2 = sa3
        sb12 = mmul(sb1,sb2)
        sa4 = mmul(ssa,sa123)
        sb3 = sa4
        sb123 = mmul(sb12,sb3)

        sb4 = mmul(ssb,sb123)
        ! Meggie out
        if ((ntest2 > 0) .and. (sa3 < sa4)) cycle

        !1.3 def mvec,c%d and c%i

        ia = a%i(sa1,sa2,sa3)
        ib = b%i(sb1,sb2,sb3)
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
        nhelp3 = dimm(b%d(0,4),sb4)

        ! sum
        nhelp41 = dimm(a%d(0,2),sa2)
        nhelp42 = dimm(a%d(0,3),sa3)
        nhelp43 = dimm(a%d(0,4),sa4)
        if ((a%d(0,6) == 2) .and. (sa2 == sa3)) then
          nhelp4 = nhelp41*(nhelp41-1)*nhelp43/2
        else if ((a%d(0,6) == 3) .and. (sa3 == sa4)) then
          nhelp4 = nhelp41*nhelp42*(nhelp42-1)/2
        else
          nhelp4 = nhelp41*nhelp42*nhelp43
        end if

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

else if (pbar == 2) then

  ! structure A(pq,rs)*B(rs,tu)=C(pq,tu)

  !2.1 define limitations - A p>q,r,s must be tested - ntest1
  !                       - A p,q,r>s must be tested - ntest2
  !                       - B r,s,t>u must be tested - ntest3

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

  if ((b%d(0,6) == 3) .or. (b%d(0,6) == 4)) then
    ntest3 = 1
  else
    ntest3 = 0
  end if

  !2.0 prepare c%d,c%i

  if ((ntest1 == 1) .and. (ntest3 == 1)) then
    nhelp1 = 4
  else if (ntest1 == 1) then
    nhelp1 = 1
  else if (ntest3 == 1) then
    nhelp1 = 3
  else
    nhelp1 = 0
  end if

  call cct3_grc0(4,nhelp1,a%d(0,1),a%d(0,2),b%d(0,3),b%d(0,4),mmul(ssa,ssb),c,posct)

  !2.2 def symm states and test the limitations

  ix = 0
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

        do sb3=1,nsym
          sb123 = mmul(sb12,sb3)

          sb4 = mmul(ssb,sb123)
          ! Meggie out
          if ((ntest3 == 1) .and. (sb3 < sb4)) cycle

          !2.3 def mvec,c%d and c%i

          ia = a%i(sa1,sa2,sa3)
          ib = b%i(sb1,sb2,sb3)
          ic = c%i(sa1,sa2,sb3)

          ! yes/no
          if ((a%d(ia,2) > 0) .and. (b%d(ib,2) > 0)) then
            nhelp1 = 1
          else
            cycle
          end if

          ! rowA
          nhelp21 = dimm(a%d(0,1),sa1)
          nhelp22 = dimm(a%d(0,2),sa2)
          if ((ntest1 == 1) .and. (sa1 == sa2)) then
            nhelp2 = nhelp21*(nhelp21-1)/2
          else
            nhelp2 = nhelp21*nhelp22
          end if

          ! colB
          nhelp31 = dimm(b%d(0,3),sb3)
          nhelp32 = dimm(b%d(0,4),sb4)
          if ((ntest3 == 1) .and. (sb3 == sb4)) then
            nhelp3 = nhelp31*(nhelp31-1)/2
          else
            nhelp3 = nhelp31*nhelp32
          end if

          ! sum
          nhelp41 = dimm(a%d(0,3),sa3)
          nhelp42 = dimm(a%d(0,4),sa4)
          if ((ntest2 == 1) .and. (sa3 == sa4)) then
            nhelp4 = nhelp41*(nhelp41-1)/2
          else
            nhelp4 = nhelp41*nhelp42
          end if

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
  end do

else if (pbar == 3) then

  ! structure A(pqr,s)*B(s,tuv)=C(pqr,tuv)
  ! not used

end if

return

end subroutine cct3_grc44c
