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

subroutine expand(wrk,wrksize,nind,exptyp,mapda,mapia,ssa,possb0,mapdb,mapib,rc)
! this routine realizes expansion
!
! A(pqrs) -> B(pqrs)
!
! nind   - # of indexex in matrix A  (Input)
! exptyp - type of expansion :  (Input)
!          1 - pq,r,s -> p,q,r,s  pq,r -> p,q,r  pq -> p,q
!          2 - p,qr,s -> p,q,r,s  p,qr -> p,q,r
!          3 - p,q,rs -> p,q,r,s
!          4 - pq,rs  -> p,q,r,s
!          5 - pq,rs  -> p,q,rs
!          6 - pq,rs  -> pq,r,s
! mapda  - direct map matrix corresponding to A  (Input)
! mapia  - inverse map matrix corresponding to A  (Input)
! ssa    - overall symmetry state  of matrix A  (Input)
! possb0 - initial position of matrix B in WRK  (Input)
! mapdb  - direct map matrix corresponding to B  (Output)
! mapib  - inverse map matrix corresponding to B  (Output)
! rc     - return (error) code  (Output)
!
! Table of expansions
!
! nind  exptyp          Operation             Implementation
! 4       0     A(p,q,r,s) -> B(p,q,r,s)     Realized in map
! 4       1     A(pq,r,s)  -> B(p,q,r,s)           Yes
! 4       2     A(p,qr,s)  -> B(p,q,r,s)           Yes
! 4       3     A(p,q,rs)  -> B(p,q,r,s)           Yes
! 4       4     A(pq,rs)   -> B(p,q,r,s)           Yes
! 4       5     A(pq,rs)   -> B(p,q,rs)            Yes
! 4       6     A(pq,rs)   -> B(pq,r,s)            Yes
!
! 3       0     A(p,q,r)   -> B(p,q,r)       Realized in map
! 3       1     A(pq,r)    -> B(p,q,r)             Yes
! 3       2     A(p,qr)    -> B(p,q,r)             Yes
!
! 2       0     A(p,q)     -> B(p,q)         Realized in map
! 2       1     A(pq)      -> B(p,q)               Yes
!
! 1       0     A(p)       -> B(p)           Realized in map

use ccsd_global, only: dimm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, nind, exptyp, mapda(0:512,6), mapdb(0:512,6), ssa, possb0, mapia(8,8,8), mapib(8,8,8), rc
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: ia, ib1, ib2, ib3, ib4, na, nhelp1, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6, posst, sa1, sa2, sa3, sa4, typa

rc = 0
na = mapda(0,5)
typa = mapda(0,6)
! general tests
if (exptyp == 0) then
  ! RC=1  : exptyp=0 (for exptyp=0, there is no sopystical expansion, NCI)
  rc = 1
  return
end if

! get mapdb,mapib

if ((nind == 4) .and. (exptyp == 5)) then
  nhelp1 = 3
else if ((nind == 4) .and. (exptyp == 6)) then
  nhelp1 = 1
else
  nhelp1 = 0
end if

call grc0(nind,nhelp1,mapda(0,1),mapda(0,2),mapda(0,3),mapda(0,4),ssa,possb0,posst,mapdb,mapib)

if (nind < 2) then
  ! RC=2 - number of indices < 2 (NCI)
  rc = 2
end if

if (nind == 2) then

  ! ********** 2 index *********

  if (exptyp == 1) then

    !2.1 expand A(pq) -> B(p,q)

    !    tests

    if (typa /= 1) then
      ! RC=3 : nind=2, exptyp=1 (typA is not 1, Stup)
      rc = 3
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)

      ib1 = mapib(sa1,1,1)
      ib2 = mapib(sa2,1,1)

      if (sa1 > sa2) then

        ! map A(p,q) -> B(p,q)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q) -> - B(q,p)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        ! dimp,dimq
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)

        call map21(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,2,1,-1)

      else
        ! sa1=sa2

        ! expand A(pq) -> B(p,q)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimp
        nhelp3 = dimm(mapda(0,1),sa1)

        call expand0(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),nhelp3)

      end if

    end do

  else
    ! RC=4 : nind=2, exptyp=@ (Improper exptyp, Stup)
    rc = 4
    return
  end if

else if (nind == 3) then

  ! ********** 3 index *********

  if (exptyp == 1) then

    !3.1 expand A(pq,r) -> B(p,q,r)

    !    tests

    if (typa /= 1) then
      ! RC=5 : nind=3, exptyp=1 (typA is not 1, Stup)
      rc = 5
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)

      ib1 = mapib(sa1,sa2,1)
      ib2 = mapib(sa2,sa1,1)

      if (sa1 > sa2) then

        ! map A(p,q,r) -> B(p,q,r)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r) -> - B(q,p,r)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        ! dimp,dimq,dimr
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)

        call map31(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,2,1,3,-1)

      else
        ! sa1=sa2

        ! expand A(pq,r) -> B(p,q,r)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimpq
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimp,dimr
        nhelp4 = dimm(mapda(0,1),sa1)
        nhelp5 = dimm(mapda(0,3),sa3)

        call expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)

      end if

    end do

  else if (exptyp == 2) then

    !3.2 expand A(p,qr) -> B(p,q,r)

    !    tests

    if (typa /= 2) then
      ! RC=6 : nind=3, exptyp=2 (typA is not 2, Stup)
      rc = 6
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)

      ib1 = mapib(sa1,sa2,1)
      ib2 = mapib(sa1,sa3,1)

      if (sa2 > sa3) then

        ! map A(p,q,r) -> B(p,q,r)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r) -> - B(q,p,r)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        ! dimp,dimq,dimr
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)

        call map31(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,1,3,2,-1)

      else
        ! sa2=sa3

        ! expand A(p,qr) -> B(p,q,r)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimp,dimq
        nhelp4 = dimm(mapda(0,1),sa1)
        nhelp5 = dimm(mapda(0,2),sa2)

        call expand3(wrk(nhelp1),wrk(nhelp2),nhelp4,nhelp5)

      end if

    end do

  else
    ! RC=7 : nind=3, exptyp=@ (Improper exptyp, Stup)
    rc = 7
    return
  end if

else if (nind == 4) then

  ! ********** 4 index *********

  if (exptyp == 1) then

    !4.1 expand A(pq,r,s) -> B(p,q,r,s)

    !    tests

    if (typa /= 1) then
      ! RC=8 : nind=4, exptyp=1 (typA is not 1, Stup)
      rc = 8
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)
      sa4 = mapda(ia,6)

      ib1 = mapib(sa1,sa2,sa3)
      ib2 = mapib(sa2,sa1,sa3)

      if (sa1 > sa2) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)

        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,3,4,-1)

      else
        ! sa1=sa2

        ! expand A(pq,r_s) -> B(p,q,r_s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimpq
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimp,dimr_s
        nhelp4 = dimm(mapda(0,1),sa1)
        nhelp5 = dimm(mapda(0,3),sa3)*dimm(mapda(0,4),sa4)

        call expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)

      end if

    end do

  else if (exptyp == 2) then

    !4.2 expand A(p,qr,s) -> B(p,q,r,s)

    !    tests

    if (typa /= 2) then
      ! RC=9 : nind=4, exptyp=2 (typA is not 2, Stup)
      rc = 9
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)
      sa4 = mapda(ia,6)

      ib1 = mapib(sa1,sa2,sa3)
      ib2 = mapib(sa1,sa3,sa2)

      if (sa2 > sa3) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r,s) -> - B(p,r,q,s)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)

        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,3,2,4,-1)

      else
        ! sa2=sa3

        ! expand A(p,qr,s) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimqr
        nhelp3 = dimm(mapda(0,2),sa1)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimp,dims,dimq
        nhelp4 = dimm(mapda(0,1),sa1)
        nhelp5 = dimm(mapda(0,4),sa4)
        nhelp5 = dimm(mapda(0,2),sa2)

        call expand2(wrk(nhelp1),wrk(nhelp2),nhelp4,nhelp3,nhelp5,nhelp6)

      end if

    end do

  else if (exptyp == 3) then

    !4.3 expand A(p,q,rs) -> B(p,q,r,s)

    !    tests

    if (typa /= 3) then
      ! RC=10: nind=4, exptyp=3 (typA is not 3, Stup)
      rc = 10
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)
      sa4 = mapda(ia,6)

      ib1 = mapib(sa1,sa2,sa3)
      ib2 = mapib(sa1,sa2,sa4)

      if (sa3 > sa4) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)

        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,2,4,3,-1)

      else
        ! sa1=sa2

        ! expand A(p_q,rs) -> B(p_q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimr,dimp_q
        nhelp4 = dimm(mapda(0,3),sa3)
        nhelp5 = dimm(mapda(0,1),sa1)*dimm(mapda(0,2),sa2)

        call expand3(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp4)

      end if

    end do

  else if (exptyp == 4) then

    !4.4 expand A(pq,rs) -> B(p,q,r,s)

    !    tests

    if (typa /= 4) then
      ! RC=11: nind=4, exptyp=4 (typA is not 4, Stup)
      rc = 11
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)
      sa4 = mapda(ia,6)

      ib1 = mapib(sa1,sa2,sa3)
      ib2 = mapib(sa2,sa1,sa3)
      ib3 = mapib(sa1,sa2,sa4)
      ib4 = mapib(sa2,sa1,sa4)

      if ((sa1 > sa2) .and. (sa3 > sa4)) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)
        ! map A(p,q,r,s) -> - B(p,q,s,r)
        ! map A(p,q,r,s) -> + B(q,p,s,r)

        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,3,4,-1)

        ! possB3
        nhelp2 = mapdb(ib3,1)
        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,2,4,3,-1)

        ! possB4
        nhelp2 = mapdb(ib4,1)
        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,4,3,1)

      else if ((sa1 == sa2) .and. (sa3 == sa4)) then

        ! expand A(pq,rs) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimp,dimq
        nhelp5 = dimm(mapda(0,1),sa1)
        nhelp6 = dimm(mapda(0,3),sa3)
        ! dimpq,dimrs
        nhelp3 = nhelp5*(nhelp5-1)/2
        nhelp4 = nhelp6*(nhelp6-1)/2

        call expand40(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6)

      else if (sa1 == sa2) then

        ! expand A(pq,r_s) -> B(p,q,r_s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! dimr,dims
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)
        ! dimpq
        nhelp4 = dimm(mapda(0,1),sa1)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! possB1
        nhelp2 = mapdb(ib1,1)
        call expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5*nhelp6,nhelp4)

        ! expand A(pq,r,s) -> - B(p,q,s,r)

        ! possB3
        nhelp2 = mapdb(ib3,1)
        call expand41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,nhelp6,nhelp4)

      else if (sa3 == sa4) then

        ! expand A(p_q,rs) -> B(p_q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! dimp,dimq
        nhelp5 = dimm(mapda(0,1),sa1)
        nhelp6 = dimm(mapda(0,2),sa2)
        ! dimrs
        nhelp4 = dimm(mapda(0,3),sa3)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! possB1
        nhelp2 = mapdb(ib1,1)
        call expand3(wrk(nhelp1),wrk(nhelp2),nhelp5*nhelp6,nhelp4)

        ! expand A(p,q,rs) -> - B(q,p,r,s)

        ! possB4
        nhelp2 = mapdb(ib2,1)
        call expand41(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,nhelp3,nhelp4)

      end if

    end do

  else if (exptyp == 5) then

    !4.5 expand A(pq,rs) -> B(p,q,rs)

    !    tests

    if (typa /= 4) then
      ! RC=12: nind=4, exptyp=5 (typA is not 4, Stup)
      rc = 12
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)
      sa4 = mapda(ia,6)

      ib1 = mapib(sa1,sa2,sa3)
      ib2 = mapib(sa2,sa1,sa3)

      if (sa3 > sa4) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)

        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)

        ! possB2
        nhelp2 = mapdb(ib2,1)
        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,3,4,-1)

      else if ((sa1 == sa2) .and. (sa3 == sa4)) then

        ! expand A(pq,rs) -> B(p,q,rs)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimp,dimr
        nhelp5 = dimm(mapda(0,1),sa1)
        nhelp6 = dimm(mapda(0,3),sa3)
        ! dimpq,dimrs
        nhelp3 = nhelp5*(nhelp5-1)/2
        nhelp4 = nhelp6*(nhelp6-1)/2

        call expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)

      else if (sa1 == sa2) then

        ! expand A(pq,r_s) -> B(p,q,r_s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! dimr,dims
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)
        ! dimpq
        nhelp4 = dimm(mapda(0,1),sa1)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! possB1
        nhelp2 = mapdb(ib1,1)
        call expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5*nhelp6,nhelp4)

      else if (sa3 == sa4) then

        ! map A(p,q,rs) -> B(p,q,rs)

        ! possA
        nhelp1 = mapda(ia,1)
        ! dimp,dimq
        nhelp5 = dimm(mapda(0,1),sa1)
        nhelp6 = dimm(mapda(0,2),sa2)
        ! dimrs
        nhelp4 = dimm(mapda(0,3),sa3)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map31(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,nhelp3,1,2,3,1)

        ! map A(p,q,rs) -> - B(q,p,rs)

        ! possB4
        nhelp2 = mapdb(ib2,1)
        call map31(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,nhelp3,2,1,3,-1)

      end if

    end do

  else if (exptyp == 6) then

    !4.6 expand A(pq,rs) -> B(pq,r,s)

    !    tests

    if (typa /= 4) then
      ! RC=13: nind=4, exptyp=6 (typA is not 4, Stup)
      rc = 13
      return
    end if

    do ia=1,na

      sa1 = mapda(ia,3)
      sa2 = mapda(ia,4)
      sa3 = mapda(ia,5)
      sa4 = mapda(ia,6)

      ib1 = mapib(sa1,sa2,sa3)
      ib3 = mapib(sa1,sa2,sa4)

      if ((sa1 > sa2) .and. (sa3 > sa4)) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map11(wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)
        ! map A(p,q,r,s) -> - B(p,q,s,r)
        ! map A(p,q,r,s) -> + B(q,p,s,r)

        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,2),sa2)
        nhelp5 = dimm(mapda(0,3),sa3)
        nhelp6 = dimm(mapda(0,4),sa4)

        ! possB3
        nhelp2 = mapdb(ib3,1)
        call map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,2,4,3,-1)

      else if ((sa1 == sa2) .and. (sa3 == sa4)) then

        ! expand A(pq,rs) -> B(pq,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! possB
        nhelp2 = mapdb(ib1,1)
        ! dimp,dimq
        nhelp5 = dimm(mapda(0,1),sa1)
        nhelp6 = dimm(mapda(0,3),sa3)
        ! dimpq,dimrs
        nhelp3 = nhelp5*(nhelp5-1)/2

        !@@? call expand4(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6)
        call expand3(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp6)

      else if (sa1 == sa2) then

        ! map A(pq,r,s) -> B(pq,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! dimp,dimr,dims
        nhelp3 = dimm(mapda(0,1),sa1)
        nhelp4 = dimm(mapda(0,3),sa3)
        nhelp5 = dimm(mapda(0,4),sa4)
        ! dimpq
        nhelp6 = nhelp3*(nhelp3-1)/2

        ! possB1
        nhelp2 = mapdb(ib1,1)
        call map31(wrk(nhelp1),wrk(nhelp2),nhelp6,nhelp4,nhelp5,1,2,3,1)

        ! map A(pq,r,s) -> - B(pq,s,r)

        ! possB3
        nhelp2 = mapdb(ib3,1)
        call map31(wrk(nhelp1),wrk(nhelp2),nhelp6,nhelp4,nhelp5,1,3,2,-1)

      else if (sa3 == sa4) then

        ! expand A(p,q,rs) -> B(p,q,r,s)

        ! possA
        nhelp1 = mapda(ia,1)
        ! dimp,dimq
        nhelp5 = dimm(mapda(0,1),sa1)
        nhelp6 = dimm(mapda(0,2),sa2)
        ! dimrs
        nhelp4 = dimm(mapda(0,3),sa3)

        ! possB1
        nhelp2 = mapdb(ib1,1)
        call expand3(wrk(nhelp1),wrk(nhelp2),nhelp5*nhelp6,nhelp4)

      end if

    end do

  else
    ! RC=14: nind=4, exptyp=@ (Improper exptyp, Stup)
    rc = 14
  end if

else
  ! RC=15: nind=@ (number of indices >4, NCI)
  rc = 15
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(mapia)

end subroutine expand
