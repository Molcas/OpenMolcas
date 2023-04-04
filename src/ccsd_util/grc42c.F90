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

subroutine grc42C(mapda,mapdb,mapdc,mapia,mapib,mapic,mvec,ssa,ssb,pbar,possc0,ix)

use ccsd_global, only: dimm, mmul, nsym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: mapda(0:512,6), mapdb(0:512,6), mapdc(0:512,6), mapia(8,8,8), mapib(8,8,8), mapic(8,8,8), mvec(4096,7), ssa, &
                     ssb, pbar, possc0, ix
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp21, nhelp22, nhelp23, nhelp3, nhelp4, nsyma2, nsyma3, ntest1, ntest2, &
                     possct, sa1, sa12, sa123, sa2, sa3, sa4, sb1, sb2

!1*

if (pbar == 1) then

  ! not possible

else if (pbar == 2) then

  ! structure A(pq,rs)*B(rs)=C(pq)
  ! implemented in grc42y

else if (pbar == 3) then

  ! structure A(pqr,s)*B(s,t)=C(pqr,t)

  !1.1 define limitations - p>q,r,s must be tested - ntest1
  !                       - p,q>r,s must be tested - ntest2

  !1.0 prepare mapdc,mapic

  call grc0(4,mapda(0,6),mapda(0,1),mapda(0,2),mapda(0,3),mapdb(0,2),mmul(ssa,ssb),possc0,possct,mapdc,mapic)

  if (mapda(0,6) == 1) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if (mapda(0,6) == 2) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

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

        !1.3 def mvec,mapdc and mapdi

        ia = mapia(sa1,sa2,sa3)
        ib = mapib(sb1,1,1)
        ic = mapic(sa1,sa2,sa3)

        ! yes/no
        if ((mapda(ia,2) <= 0) .or. (mapdb(ib,2) <= 0)) cycle
        nhelp1 = 1

        ! rowA
        nhelp21 = dimm(mapda(0,1),sa1)
        nhelp22 = dimm(mapda(0,2),sa2)
        nhelp23 = dimm(mapda(0,3),sa3)
        if ((ntest1 == 1) .and. (sa1 == sa2)) then
          nhelp2 = nhelp21*(nhelp21-1)*nhelp23/2
        else if ((ntest2 == 1) .and. (sa2 == sa3)) then
          nhelp2 = nhelp21*nhelp22*(nhelp22-1)/2
        else
          nhelp2 = nhelp21*nhelp22*nhelp23
        end if

        ! colB
        nhelp3 = dimm(mapdb(0,2),sb2)

        ! sum
        nhelp4 = dimm(mapda(0,4),sa4)

        mvec(ix,1) = nhelp1
        mvec(ix,2) = mapda(ia,1)
        mvec(ix,3) = mapdb(ib,1)
        mvec(ix,4) = mapdc(ic,1)
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

end subroutine grc42C
