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

subroutine cct3_expand(wrk,wrksize,nind,exptyp,a,ssa,b,rc)
! this routine realizes expansion
!
! A(pqrs) -> B(pqrs)
!
! nind   - # of indices in matrix A  (Input)
! exptyp - type of expansion :  (Input)
!          1 - pq,r,s -> p,q,r,s  pq,r -> p,q,r  pq -> p,q
!          2 - p,qr,s -> p,q,r,s  p,qr -> p,q,r
!          3 - p,q,rs -> p,q,r,s
!          4 - pq,rs  -> p,q,r,s
!          5 - pq,rs  -> p,q,rs
!          6 - pq,rs  -> pq,r,s
! a     - A  (Input)
! ssa   - overall symmetry state of matrix A  (Input)
! b     - B  (Input/Output)
! rc    - return (error) code  (Output)
!
! Table of expansions
!
! nind  exptyp          Operation             Implementation
!
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

use CCT3_global, only: dimm, Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, nind, exptyp, ssa
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: a
type(Map_Type), intent(inout) :: b
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: ia, ib1, ib2, ib3, ib4, na, nhelp1, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6, post, sa1, sa2, sa3, sa4, typa

rc = 0
na = a%d(0,5)
typa = a%d(0,6)

! general tests

if (exptyp == 0) then
  ! RC=1  : exptyp=0 (for exptyp=0, there is no sophistical expansion, NCI)
  rc = 1
  return
end if

! get b%d,b%i

if ((nind == 4) .and. (exptyp == 5)) then
  nhelp1 = 3
else if ((nind == 4) .and. (exptyp == 6)) then
  nhelp1 = 1
else
  nhelp1 = 0
end if

call cct3_grc0(nind,nhelp1,a%d(0,1),a%d(0,2),a%d(0,3),a%d(0,4),ssa,b,post)

if (nind < 2) then
  ! RC=2 - number of indices < 2 (NCI)
  rc = 2
end if

if (nind == 2) then

  ! ********** 2 index *********

  if (exptyp == 1) then

    !2.1 expand A(pq) -> B(p,q)

    ! tests

    if (typa /= 1) then
      ! RC=3 : nind=2, exptyp=1 (typA is not 1, Stup)
      rc = 3
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)

      ib1 = b%i(sa1,1,1)
      ib2 = b%i(sa2,1,1)

      if (sa1 > sa2) then

        ! map A(p,q) -> B(p,q)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q) -> - B(q,p)

        ! posB2
        nhelp2 = b%d(ib2,1)
        ! dimp,dimq
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)

        call cct3_map21(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,2,1,-1)

      else
        ! sa1=sa2

        ! expand A(pq) -> B(p,q)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimp
        nhelp3 = dimm(a%d(0,1),sa1)

        call cct3_expand0(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),nhelp3)

      end if

    end do

  else
    ! RC=4 : nind=2, exptyp=@ (Inproper exptyp, Stup)
    rc = 4
    return
  end if

else if (nind == 3) then

  ! ********** 3 index *********

  if (exptyp == 1) then

    !3.1 expand A(pq,r) -> B(p,q,r)

    ! tests

    if (typa /= 1) then
      ! RC=5 : nind=3, exptyp=1 (typA is not 1, Stup)
      rc = 5
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)

      ib1 = b%i(sa1,sa2,1)
      ib2 = b%i(sa2,sa1,1)

      if (sa1 > sa2) then

        ! map A(p,q,r) -> B(p,q,r)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r) -> - B(q,p,r)

        ! posB2
        nhelp2 = b%d(ib2,1)
        ! dimp,dimq,dimr
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)

        call cct3_map31(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,2,1,3,-1)

      else
        ! sa1=sa2

        ! expand A(pq,r) -> B(p,q,r)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimpq
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimp,dimr
        nhelp4 = dimm(a%d(0,1),sa1)
        nhelp5 = dimm(a%d(0,3),sa3)

        call cct3_expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)

      end if

    end do

  else if (exptyp == 2) then

    !3.2 expand A(p,qr) -> B(p,q,r)

    ! tests

    if (typa /= 2) then
      ! RC=6 : nind=3, exptyp=2 (typA is not 2, Stup)
      rc = 6
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)

      ib1 = b%i(sa1,sa2,1)
      ib2 = b%i(sa1,sa3,1)

      if (sa2 > sa3) then

        ! map A(p,q,r) -> B(p,q,r)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r) -> - B(q,p,r)

        ! posB2
        nhelp2 = b%d(ib2,1)
        ! dimp,dimq,dimr
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)

        call cct3_map31(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,1,3,2,-1)

      else
        ! sa2=sa3

        ! expand A(p,qr) -> B(p,q,r)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimqr
        nhelp3 = dimm(a%d(0,2),sa1)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimp,dimq
        nhelp4 = dimm(a%d(0,1),sa1)
        nhelp5 = dimm(a%d(0,2),sa2)

        call cct3_expand3(wrk(nhelp1),wrk(nhelp2),nhelp4,nhelp3,nhelp5)

      end if

    end do

  else
    ! RC=7 : nind=3, exptyp=@ (Inproper exptyp, Stup)
    rc = 7
    return
  end if

else if (nind == 4) then

  ! ********** 4 index *********

  if (exptyp == 1) then

    !4.1 expand A(pq,r,s) -> B(p,q,r,s)

    ! tests

    if (typa /= 1) then
      ! RC=8 : nind=4, exptyp=1 (typA is not 1, Stup)
      rc = 8
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)
      sa4 = a%d(ia,6)

      ib1 = b%i(sa1,sa2,sa3)
      ib2 = b%i(sa2,sa1,sa3)

      if (sa1 > sa2) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)

        ! posB2
        nhelp2 = b%d(ib2,1)
        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)

        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,3,4,-1)

      else
        ! sa1=sa2

        ! expand A(pq,r_s) -> B(p,q,r_s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimpq
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimp,dimr_s
        nhelp4 = dimm(a%d(0,1),sa1)
        nhelp5 = dimm(a%d(0,3),sa3)*dimm(a%d(0,4),sa4)

        call cct3_expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)

      end if

    end do

  else if (exptyp == 2) then

    !4.2 expand A(p,qr,s) -> B(p,q,r,s)

    ! tests

    if (typa /= 2) then
      ! RC=9 : nind=4, exptyp=2 (typA is not 2, Stup)
      rc = 9
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)
      sa4 = a%d(ia,6)

      ib1 = b%i(sa1,sa2,sa3)
      ib2 = b%i(sa1,sa3,sa2)

      if (sa2 > sa3) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r,s) -> - B(p,r,q,s)

        ! posB2
        nhelp2 = b%d(ib2,1)
        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)

        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,3,2,4,-1)

      else
        ! sa2=sa3

        ! expand A(p,qr,s) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimqr
        nhelp3 = dimm(a%d(0,2),sa1)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimp,dims,dimq
        nhelp4 = dimm(a%d(0,1),sa1)
        nhelp5 = dimm(a%d(0,4),sa4)
        nhelp5 = dimm(a%d(0,2),sa2)

        call cct3_expand2(wrk(nhelp1),wrk(nhelp2),nhelp4,nhelp3,nhelp5,nhelp6)

      end if

    end do

  else if (exptyp == 3) then

    !4.3 expand A(p,q,rs) -> B(p,q,r,s)

    ! tests

    if (typa /= 3) then
      ! RC=10: nind=4, exptyp=3 (typA is not 3, Stup)
      rc = 10
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)
      sa4 = a%d(ia,6)

      ib1 = b%i(sa1,sa2,sa3)
      ib2 = b%i(sa1,sa2,sa4)

      if (sa3 > sa4) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)

        ! posB2
        nhelp2 = b%d(ib2,1)
        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)

        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,2,4,3,-1)

      else
        ! sa1=sa2

        ! expand A(p_q,rs) -> B(p_q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimrs
        nhelp3 = dimm(a%d(0,3),sa3)
        nhelp3 = nhelp3*(nhelp3-1)/2
        ! dimr,dimp_q
        nhelp4 = dimm(a%d(0,3),sa3)
        nhelp5 = dimm(a%d(0,1),sa1)*dimm(a%d(0,2),sa2)

        call cct3_expand3(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp3,nhelp4)

      end if

    end do

  else if (exptyp == 4) then

    !4.4 expand A(pq,rs) -> B(p,q,r,s)

    ! tests

    if (typa /= 4) then
      ! RC=11: nind=4, exptyp=4 (typA is not 4, Stup)
      rc = 11
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)
      sa4 = a%d(ia,6)

      ib1 = b%i(sa1,sa2,sa3)
      ib2 = b%i(sa2,sa1,sa3)
      ib3 = b%i(sa1,sa2,sa4)
      ib4 = b%i(sa2,sa1,sa4)

      if ((sa1 > sa2) .and. (sa3 > sa4)) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)
        ! map A(p,q,r,s) -> - B(p,q,s,r)
        ! map A(p,q,r,s) -> + B(q,p,s,r)

        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)

        ! posB2
        nhelp2 = b%d(ib2,1)
        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,3,4,-1)

        ! posB3
        nhelp2 = b%d(ib3,1)
        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,2,4,3,-1)

        ! posB4
        nhelp2 = b%d(ib4,1)
        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,4,3,1)

      else if ((sa1 == sa2) .and. (sa3 == sa4)) then

        ! expand A(pq,rs) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimp,dimq
        nhelp5 = dimm(a%d(0,1),sa1)
        nhelp6 = dimm(a%d(0,3),sa3)
        ! dimpq,dimrs
        nhelp3 = nhelp5*(nhelp5-1)/2
        nhelp4 = nhelp6*(nhelp6-1)/2

        call cct3_expand40(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6)

      else if (sa1 == sa2) then

        ! expand A(pq,r_s) -> B(p,q,r_s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! dimr,dims
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)
        ! dimpq
        nhelp4 = dimm(a%d(0,1),sa1)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5*nhelp6,nhelp4)

        ! expand A(pq,r,s) -> - B(p,q,s,r)

        ! posB3
        nhelp2 = b%d(ib3,1)
        call cct3_expand41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,nhelp6,nhelp4)

      else if (sa3 == sa4) then

        ! expand A(p_q,rs) -> B(p_q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! dimp,dimq
        nhelp5 = dimm(a%d(0,1),sa1)
        nhelp6 = dimm(a%d(0,2),sa2)
        ! dimrs
        nhelp4 = dimm(a%d(0,3),sa3)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_expand3(wrk(nhelp1),wrk(nhelp2),nhelp5*nhelp6,nhelp3,nhelp4)

        ! expand A(p,q,rs) -> - B(q,p,r,s)

        ! posB4
        nhelp2 = b%d(ib2,1)
        call cct3_expand41(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,nhelp3,nhelp4)

      end if

    end do

  else if (exptyp == 5) then

    !4.5 expand A(pq,rs) -> B(p,q,rs)

    ! tests

    if (typa /= 4) then
      ! RC=12: nind=4, exptyp=5 (typA is not 4, Stup)
      rc = 12
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)
      sa4 = a%d(ia,6)

      ib1 = b%i(sa1,sa2,sa3)
      ib2 = b%i(sa2,sa1,sa3)

      if (sa3 > sa4) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)

        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)

        ! posB2
        nhelp2 = b%d(ib2,1)
        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,2,1,3,4,-1)

      else if ((sa1 == sa2) .and. (sa3 == sa4)) then

        ! expand A(pq,rs) -> B(p,q,rs)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimp,dimr
        nhelp5 = dimm(a%d(0,1),sa1)
        nhelp6 = dimm(a%d(0,3),sa3)
        ! dimpq,dimrs
        nhelp3 = nhelp5*(nhelp5-1)/2
        nhelp4 = nhelp6*(nhelp6-1)/2

        call cct3_expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)

      else if (sa1 == sa2) then

        ! expand A(pq,r_s) -> B(p,q,r_s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! dimr,dims
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)
        ! dimpq
        nhelp4 = dimm(a%d(0,1),sa1)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_expand1(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5*nhelp6,nhelp4)

      else if (sa3 == sa4) then

        ! map A(p,q,rs) -> B(p,q,rs)

        ! posA
        nhelp1 = a%d(ia,1)
        ! dimp,dimq
        nhelp5 = dimm(a%d(0,1),sa1)
        nhelp6 = dimm(a%d(0,2),sa2)
        ! dimrs
        nhelp4 = dimm(a%d(0,3),sa3)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map31(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,nhelp3,1,2,3,1)

        ! map A(p,q,rs) -> - B(q,p,rs)

        ! posB4
        nhelp2 = b%d(ib2,1)
        call cct3_map31(wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,nhelp3,2,1,3,-1)

      end if

    end do

  else if (exptyp == 6) then

    !4.6 expand A(pq,rs) -> B(pq,r,s)

    ! tests

    if (typa /= 4) then
      ! RC=13: nind=4, exptyp=6 (typA is not 4, Stup)
      rc = 13
      return
    end if

    do ia=1,na

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)
      sa4 = a%d(ia,6)

      ib1 = b%i(sa1,sa2,sa3)
      ib3 = b%i(sa1,sa2,sa4)

      if ((sa1 > sa2) .and. (sa3 > sa4)) then

        ! map A(p,q,r,s) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map11(wrk(nhelp1),wrk(nhelp2),a%d(ia,2),1)

        ! map A(p,q,r,s) -> - B(q,p,r,s)
        ! map A(p,q,r,s) -> - B(p,q,s,r)
        ! map A(p,q,r,s) -> + B(q,p,s,r)

        ! dimp,dimq,dimr,dims
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,2),sa2)
        nhelp5 = dimm(a%d(0,3),sa3)
        nhelp6 = dimm(a%d(0,4),sa4)

        ! posB3
        nhelp2 = b%d(ib3,1)
        call cct3_map41(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6,1,2,4,3,-1)

      else if ((sa1 == sa2) .and. (sa3 == sa4)) then

        ! expand A(pq,rs) -> B(pq,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! posB
        nhelp2 = b%d(ib1,1)
        ! dimp,dimq
        nhelp5 = dimm(a%d(0,1),sa1)
        nhelp6 = dimm(a%d(0,3),sa3)
        ! dimpq,dimrs
        nhelp3 = nhelp5*(nhelp5-1)/2
        nhelp4 = nhelp6*(nhelp6-1)/2

        !@@? call cct3_expand4 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6)
        call cct3_expand3(wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp6)

      else if (sa1 == sa2) then

        ! map A(pq,r,s) -> B(pq,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! dimp,dimr,dims
        nhelp3 = dimm(a%d(0,1),sa1)
        nhelp4 = dimm(a%d(0,3),sa3)
        nhelp5 = dimm(a%d(0,4),sa4)
        ! dimpq
        nhelp6 = nhelp3*(nhelp3-1)/2

        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_map31(wrk(nhelp1),wrk(nhelp2),nhelp6,nhelp4,nhelp5,1,2,3,1)

        ! map A(pq,r,s) -> - B(pq,s,r)

        ! posB3
        nhelp2 = b%d(ib3,1)
        call cct3_map31(wrk(nhelp1),wrk(nhelp2),nhelp6,nhelp4,nhelp5,1,3,2,-1)

      else if (sa3 == sa4) then

        ! expand A(p,q,rs) -> B(p,q,r,s)

        ! posA
        nhelp1 = a%d(ia,1)
        ! dimp,dimq
        nhelp5 = dimm(a%d(0,1),sa1)
        nhelp6 = dimm(a%d(0,2),sa2)
        ! dimrs
        nhelp4 = dimm(a%d(0,3),sa3)
        nhelp3 = nhelp4*(nhelp4-1)/2

        ! posB1
        nhelp2 = b%d(ib1,1)
        call cct3_expand3(wrk(nhelp1),wrk(nhelp2),nhelp5*nhelp6,nhelp3,nhelp4)

      end if

    end do

  else
    ! RC=14: nind=4, exptyp=@ (Inproper exptyp, Stup)
    rc = 14
  end if

else
  ! RC=15- nind=@ (number of indices >4, NCI)
  rc = 15
end if

return

end subroutine cct3_expand
