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

subroutine cct3_grc34C(mapda,mapdb,mapdc,mapia,mapib,mapic,mvec,ssa,ssb,pbar,possc0,ix)

#include "t31.fh"
integer mapda(0:512,1:6)
integer mapdb(0:512,1:6)
integer mapdc(0:512,1:6)
integer mapia(1:8,1:8,1:8)
integer mapib(1:8,1:8,1:8)
integer mapic(1:8,1:8,1:8)
integer mvec(1:4096,1:7)
integer pbar, possc0
integer ssa, ssb
! help variables
integer nhelp1, nhelp2, nhelp3, nhelp4
integer nhelp31, nhelp32
integer nhelp41, nhelp42
integer ntest1, ntest2
integer sa1, sa2, sa3, sb1, sb2, sb3, sb4, sa12, sb12, sb123
integer ia, ib, ic, ix
integer possct

!1*

if (pbar == 1) then

  ! structure A(p,qr)*B(qr,st)=C(p,st)

  !1.1  define limitations - q>r,s,t must be tested - ntest1
  !                        - q,r,s>t must be tested - ntest2

  if ((mapdb(0,6) == 1) .or. (mapdb(0,6) == 4)) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if ((mapdb(0,6) == 3) .or. (mapdb(0,6) == 4)) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

  !1.0 prepare mapdc,mapic

  if (ntest2 == 1) then
    nhelp1 = 2
  else
    nhelp1 = 0
  end if

  call cct3_grc0(3,nhelp1,mapda(0,1),mapdb(0,3),mapdb(0,4),0,mmul(ssa,ssb),possc0,possct,mapdc,mapic)

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

        !1.3 def mvec,mapdc and mapdi

        ia = mapia(sa1,sa2,1)
        ib = mapib(sb1,sb2,sb3)
        ic = mapic(sa1,sb3,1)

        ! yes/no
        if ((mapda(ia,2) > 0) .and. (mapdb(ib,2) > 0)) then
          nhelp1 = 1
        else
          cycle
        end if

        ! rowA
        nhelp2 = dimm(mapda(0,1),sa1)

        ! colB
        nhelp31 = dimm(mapdb(0,3),sb3)
        nhelp32 = dimm(mapdb(0,4),sb4)
        if ((ntest2 == 1) .and. (sb3 == sb4)) then
          nhelp3 = nhelp31*(nhelp31-1)/2
        else
          nhelp3 = nhelp31*nhelp32
        end if

        ! sum
        nhelp41 = dimm(mapda(0,2),sa2)
        nhelp42 = dimm(mapda(0,3),sa3)
        if ((ntest1 == 1) .and. (sa2 == sa3)) then
          nhelp4 = nhelp41*(nhelp41-1)/2
        else
          nhelp4 = nhelp41*nhelp42
        end if

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

else if (pbar == 3) then

  ! structure A(pq,r)*B(r,stu)=C(pq,stu)
  ! not used

end if
ix = ix-1

return

end subroutine cct3_grc34C
