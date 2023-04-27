!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Pavel Neogrady                                   *
!***********************************************************************

subroutine multstack(wrk,wrksize,a,b,c,ssa,ssb,bsize)
! This is a special routine for multiplying of the:
! C(ij,Bp) = A(ij,cd) . B(cd,Bp)
! where Bp is a limited (partial) summation over
! b index (#b - bsize), namely those, which are stacked
!
! This routine is used only in stacking in sumoverab
! process and is a modification of grc42y routine.
! Type of index Bp is registered as for standard b, but
! all lengths of blocks are calculated with
! bsize, instead of dimm(typb,symb) and the symmetry
! of b is ignored, B and C are treated as 2index
!
! P.N. 17.02.06

use ccsd_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, ssa, ssb, bsize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp) :: ia, ib, ix, iy, mvec(4096,7), nhelp1, nhelp2, nhelp21, nhelp22, nhelp3, nhelp4, nhelp41, nhelp42, ntest1, &
                     ntest2, posct, sa1, sa134, sa2, sa3, sa34, sa4, sb1, sb2

!1*

! structure A(pq,rs)*B(rs,t)=C(pq,t)
! structure A(pq,rs)*B(rs)=YC(pq)

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

call grc0stack(bsize,ntest1,a%d(0,1),a%d(0,2),b%d(0,3),0,mmul(ssa,ssb),posct,c)

!1.2 def symm states and test the limitations

ix = 1
do sb1=1,nsym
  sa3 = sb1

  sb2 = mmul(ssb,sb1)
  sa4 = sb2
  sa34 = mmul(sa3,sa4)
  ! Meggie out
  if ((ntest2 == 1) .and. (sb1 < sb2)) cycle

  do sa1=1,nsym
    sa134 = mmul(sa1,sa34)

    sa2 = mmul(ssa,sa134)
    ! Meggie out
    if ((ntest1 == 1) .and. (sa1 < sa2)) cycle

    !1.3 def mvec, c%d and c%i

    ia = a%i(sa1,sa2,sa3)
    ib = b%i(sb1,1,1)
    iy = c%i(sa1,1,1)

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

    ! sum
    nhelp41 = dimm(a%d(0,3),sa3)
    nhelp42 = dimm(a%d(0,4),sa4)
    if ((ntest2 == 1) .and. (sa3 == sa4)) then
      nhelp4 = nhelp41*(nhelp41-1)/2
    else
      nhelp4 = nhelp41*nhelp42
    end if

    ! colBp
    nhelp3 = bsize

    mvec(ix,1) = nhelp1
    mvec(ix,2) = a%d(ia,1)
    mvec(ix,3) = b%d(ib,1)
    mvec(ix,4) = c%d(iy,1)
    mvec(ix,5) = nhelp2
    mvec(ix,6) = nhelp4
    mvec(ix,7) = nhelp3

    ix = ix+1

  end do
end do
ix = ix-1

!* multiplying

call multc0(wrk,wrksize,mvec,ix,c,1)

return

end subroutine multstack
