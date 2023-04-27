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

subroutine add(wrk,wrksize,ninda,nindb,nindext,typext,u,v,ssu,ssv,factor,a,ssa,b,ssb,rc)
! this routine does:
! B(indb) = B(indb) + factor * A(inda)
!
! ninda   - # of indices in matrix A (1-4)
! nindb   - # of indices in matrix B (1-4)
! nindext - # of external (frozen, fixed) indices (now:0-2)
! typext  - characterize external indices as follows:
!           0 - no frozen index
!           1 - frozen index p
!           2 - frozen index q
!           3 - frozen index r
!           4 - frozen index s
!           5 - frozen indices p,q
!           6 - frozen indices r,s
! u       - value of first external index (if any, else 0)
! v       - value of second external index (if any, else 0)
! ssu     - symmetry of u (if any, else 1)
! ssv     - symmetry of v (if any, else 1)
! factor  - multiplicative factor (see def)
! a       - map type corresponding to A
! ssa     - overall spin state of matrix A
! b       - map type corresponding to B
! ssb     - overall spin state of matrix B
! rc      - return (error) code
!
! Table of present implementations:
!
! nindB  nindxet  typext  typB  =>   Implemented
! >4                                   No
!
! 4       0        0    0-4            Yes
! 4       1       1-4   0,4            Yes
! 4       1       1-4   2,3            No
! 4       2        5    0,4            Yes
! 4       2        5    2,3            No
! 4       2       6-n                  No
! 4       3                            No
!
! 3       0       0     0-2            Yes
! 3       1       1-3    0             Yes
! 3       1       1-3   1,2            No
! 3       2                            No
!
! 2       0       0     0,1            Yes
! 2       1       1-2    0             Yes
! 2       1       1-2    1             No
! 2       2                            No
!
! 1                                    No
!
! !N.B. oprav co je oznacene c@!

use ccsd_global, only: dimm, Map_Type, mmul
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, ninda, nindb, nindext, typext, u, v, ssu, ssv, ssa, ssb
real(kind=wp), intent(inout) :: wrk(wrksize)
real(kind=wp), intent(in) :: factor
type(Map_Type), intent(in) :: a, b
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: ia, ib, ibm, nhelp1, nhelp10, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6, nhelp7, nhelp8, nhelp9, p, pq, q, sa1, &
                     sa2, sa3, ssp, ssq, typa, typb
real(kind=wp) :: fact

! To fix compiler warnings
p = 0
q = 0
! general tests

nhelp1 = nindA+nindext

if (nhelp1 /= nindb) then
  ! RC=1  : incompatible (nindA, nindB and nindext, Stup)
  rc = 1
  return
end if

nhelp1 = mmul(ssu,ssv)
nhelp1 = mmul(ssa,nhelp1)
if (nhelp1 /= ssb) then
  ! RC=2  : incompatible (ssa, ssb ,ssu and ssv, Stup)
  rc = 2
  return
end if

typa = a%d(0,6)
typb = b%d(0,6)
fact = factor

if ((typb >= 1) .and. (typb <= 3)) then
  ! RC=3 : typB is 1,2 or 3 (NCI)
  rc = 3
  return
end if

if ((nindext == 2) .and. (typb == 4)) then

  ! def p,q,ssp,ssq,fact(new)

  if (ssu > ssv) then
    ! ssu>ssv
    p = u
    q = v
    ssp = ssu
    ssq = ssv
    fact = factor
  else if (ssu == ssv) then
    ! ssu=ssv
    if (u >= v) then
      p = u
      q = v
      ssp = ssu
      ssq = ssv
      fact = factor
    else
      p = v
      q = u
      ssp = ssv
      ssq = ssu
      fact = -factor
    end if
  else
    ! ssu<ssv
    p = v
    q = u
    ssp = ssv
    ssq = ssu
    fact = -factor
  end if

end if

if (nindb == 4) then

  ! **********  -> B(pqrs) **********

  if (nindext == 0) then

    !400  case B(pqrs) <-- A(pqrs) or B(p,q,r,s) <-- A(p,q,r,s)

    !     tests

    if (typa /= typb) then
      ! RC=4 : nindB=4, nindext=0 (TypA incompatible with TypB ,Stup)
      rc = 4
      return
    end if

    do ia=1,a%d(0,5)

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)

      ib = b%i(sa1,sa2,sa3)

      ! def length
      nhelp1 = a%d(ia,2)
      if (nhelp1 == 0) cycle

      ! def posA,posB
      nhelp2 = a%d(ia,1)
      nhelp3 = b%d(ib,1)

      call add10(wrk(nhelp2),wrk(nhelp3),nhelp1,fact)

    end do

  else if (nindext == 1) then

    if (typext == 1) then

      if (typb == 0) then

        !4110 case B(p,q,r,s) <-- A(q,r,s)

        !     tests

        if (typa /= 0) then
          ! RC=5 : nindB=4, nindeext=1, typext=1, typB=0 (typA is not 0, Stup)
          rc = 5
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(ssu,sa1,sa2)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! def fictive dimensions
          nhelp8 = nhelp5*nhelp6*nhelp7

          call add21(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp8,fact)

        end do

      else if (typb == 4) then

        !4114 case B(pq,rs) <-- A(q,rs)

        !     tests

        if (typa /= 2) then
          ! RC=6  : nindB=4, nindeext=1, typext=1, typB=4 (typA is not 2, Stup)
          rc = 6
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(ssu,sa1,sa2)
          ibm = b%i(sa1,ssu,sa2)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA
          nhelp2 = a%d(ia,1)

          if (ssu > sa1) then

            ! def dimp,dimq,dimr,dims
            nhelp4 = dimm(b%d(0,1),b%d(ib,3))
            nhelp5 = dimm(b%d(0,2),b%d(ib,4))
            nhelp6 = dimm(b%d(0,3),b%d(ib,5))
            nhelp7 = dimm(b%d(0,4),b%d(ib,6))

            ! def fictive dimensions
            if (sa2 == sa3) then
              nhelp8 = nhelp6*(nhelp6-1)/2
            else
              nhelp8 = nhelp6*nhelp7
            end if

            ! def posB
            nhelp3 = b%d(ib,1)
            ! def fictive dimensions
            nhelp9 = nhelp5*nhelp8
            call add21(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,fact)

          else if (ssu == sa1) then

            ! def dimp,dimq,dimr,dims
            nhelp4 = dimm(b%d(0,1),b%d(ib,3))
            nhelp5 = dimm(b%d(0,2),b%d(ib,4))
            nhelp6 = dimm(b%d(0,3),b%d(ib,5))
            nhelp7 = dimm(b%d(0,4),b%d(ib,6))

            ! def fictive dimensions
            if (sa2 == sa3) then
              nhelp8 = nhelp6*(nhelp6-1)/2
            else
              nhelp8 = nhelp6*nhelp7
            end if

            ! def posB
            nhelp3 = b%d(ib,1)
            ! def fictive dimensions
            nhelp9 = nhelp4*(nhelp4-1)/2
            call add41(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,nhelp8,fact)

          else
            ! ssu<sa1  B(qp,rs) <-- -A_p (q,rs)

            ! def dimp,dimq,dimr,dims
            nhelp4 = dimm(b%d(0,1),b%d(ibm,3))
            nhelp5 = dimm(b%d(0,2),b%d(ibm,4))
            nhelp6 = dimm(b%d(0,3),b%d(ibm,5))
            nhelp7 = dimm(b%d(0,4),b%d(ibm,6))

            ! def fictive dimensions
            if (sa2 == sa3) then
              nhelp8 = nhelp6*(nhelp6-1)/2
            else
              nhelp8 = nhelp6*nhelp7
            end if

            ! def posB-
            nhelp3 = b%d(ibm,1)
            call add32(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,nhelp8,-fact)

          end if

        end do

      else
        ! RC=7 : nindB=4, nindext=1, typext=1, (typA is not 0 or 4, (NCI))
        rc = 7
        return
      end if

    else if (typext == 2) then

      if (typb == 0) then

        !4120 case B(p,q,r,s) <-- A(p,r,s)

        !     tests

        if (typa /= 0) then
          ! RC=8 : nindB=4, nindeext=1, typext=2, typB=0 (typA is not 0, Stup)
          rc = 8
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(sa1,ssu,sa2)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! def fictive dimensions
          nhelp8 = nhelp6*nhelp7

          call add32(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,nhelp8,fact)

        end do

      else if (typb == 4) then

        !4124 case B(pq,rs) <-- A(p,rs)

        !     tests

        if (typa /= 2) then
          ! RC=9 : nindB=4, nindeext=1, typext=2, typB=4 (typA is not 2, Stup)
          rc = 9
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(sa1,ssu,sa2)
          ibm = b%i(ssu,sa1,sa2)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA
          nhelp2 = a%d(ia,1)

          if (sa1 > ssu) then

            ! def dimp,dimq,dimr,dims
            nhelp4 = dimm(b%d(0,1),b%d(ib,3))
            nhelp5 = dimm(b%d(0,2),b%d(ib,4))
            nhelp6 = dimm(b%d(0,3),b%d(ib,5))
            nhelp7 = dimm(b%d(0,4),b%d(ib,6))

            ! def fictive dimensions
            if (sa2 == sa3) then
              nhelp8 = nhelp6*(nhelp6-1)/2
            else
              nhelp8 = nhelp6*nhelp7
            end if

            ! def posB
            nhelp3 = b%d(ib,1)
            call add32(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,nhelp8,fact)

          else if (sa1 == ssu) then

            ! def dimp,dimq,dimr,dims
            nhelp4 = dimm(b%d(0,1),b%d(ib,3))
            nhelp5 = dimm(b%d(0,2),b%d(ib,4))
            nhelp6 = dimm(b%d(0,3),b%d(ib,5))
            nhelp7 = dimm(b%d(0,4),b%d(ib,6))

            ! def fictive dimensions
            if (sa2 == sa3) then
              nhelp8 = nhelp6*(nhelp6-1)/2
            else
              nhelp8 = nhelp6*nhelp7
            end if

            ! def posB
            nhelp3 = b%d(ib,1)
            ! def fictive dimensions
            nhelp9 = nhelp4*(nhelp4-1)/2
            call add42(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,nhelp8,fact)

          else
            ! sa1<ssu  B(qp,rs) <-- -A_q (p,rs)

            ! def dimp,dimq,dimr,dims
            nhelp4 = dimm(b%d(0,1),b%d(ibm,3))
            nhelp5 = dimm(b%d(0,2),b%d(ibm,4))
            nhelp6 = dimm(b%d(0,3),b%d(ibm,5))
            nhelp7 = dimm(b%d(0,4),b%d(ibm,6))

            ! def fictive dimensions
            if (sa2 == sa3) then
              nhelp8 = nhelp6*(nhelp6-1)/2
            else
              nhelp8 = nhelp6*nhelp7
            end if

            ! def posB-
            nhelp3 = b%d(ibm,1)
            ! def fictive index
            nhelp9 = nhelp8*nhelp5
            call add21(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp9,-fact)

          end if

        end do

      else
        ! RC=10: nindB=4, nindext=1, typext=2, (typA is not 0 or 4, (NCI))
        rc = 10
        return
      end if

    else if (typext == 3) then

      if (typb == 0) then

        !4130 case B(p,q,r,s) <-- A(p,q,s)

        !     tests

        if (typa /= 0) then
          ! RC=11: nindB=4, nindeext=1, typext=3, typB=0 (typA is not 0, Stup)
          rc = 11
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(sa1,sa2,ssu)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! def fictive dimensions
          nhelp8 = nhelp4*nhelp5

          call add32(wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp6,nhelp7,fact)

        end do

      else if (typb == 4) then

        !4134 case B(pq,rs) <-- A(pq,s)
        !@!   oprav to tak ako v typext 1 a 2

        !     tests

        if (typa /= 1) then
          ! RC=12: nindB=4, nindeext=1, typext=3, typB=4 (typA is not 1, Stup)
          rc = 12
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(sa1,sa2,ssu)
          ibm = b%i(sa1,sa2,sa3)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA
          nhelp2 = a%d(ia,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! def fictive dimensions
          if (sa1 == sa2) then
            nhelp8 = nhelp4*(nhelp4-1)/2
          else
            nhelp8 = nhelp4*nhelp5
          end if

          if (ssu > sa3) then

            ! def posB
            nhelp3 = b%d(ib,1)
            call add32(wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp6,nhelp7,fact)

          else if (ssu == sa3) then

            ! def posB
            nhelp3 = b%d(ib,1)
            ! def fictive dimensions
            nhelp9 = nhelp6*(nhelp6-1)/2
            call add43(wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp9,nhelp6,fact)

          else
            ! ssu<sa3  B(pq,sr) <-- -A_r (pq,s)

            ! def posB-
            nhelp3 = b%d(ibm,1)
            ! def fictive dimension
            nhelp9 = nhelp8*nhelp7
            call add22(wrk(nhelp2),wrk(nhelp3),u,nhelp9,nhelp6,-fact)

          end if

        end do

      else
        ! RC=13: nindB=4, nindext=1, typext=3, (typA is not 0 or 4, (NCI))
        rc = 13
        return
      end if

    else if (typext == 4) then

      if (typb == 0) then

        !4140 case B(p,q,r,s) <-- A(p,q,r)

        !     tests

        if (typa /= 0) then
          ! RC=14: nindB=4, nindeext=1, typext=4, typB=0 (typA is not 0, Stup)
          rc = 14
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(sa1,sa2,sa3)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! def fictive dimensions
          nhelp8 = nhelp4*nhelp5*nhelp6

          call add22(wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp7,fact)

        end do

      else if (typb == 4) then

        !4144 case B(pq,rs) <-- A(pq,r)
        !@!   oprav to tak ako v typext 1 a 2

        !     tests

        if (typa /= 1) then
          ! RC=15: nindB=4, nindeext=1, typext=4, typB=4 (typA is not 1, Stup)
          rc = 15
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)
          sa3 = a%d(ia,5)

          ib = b%i(sa1,sa2,sa3)
          ibm = b%i(sa1,sa2,ssu)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA
          nhelp2 = a%d(ia,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! def fictive dimensions
          if (sa1 == sa2) then
            nhelp8 = nhelp4*(nhelp4-1)/2
          else
            nhelp8 = nhelp4*nhelp5
          end if

          if (sa3 > ssu) then

            ! def posB
            nhelp3 = b%d(ib,1)
            ! def fictive dimension
            nhelp9 = nhelp8*nhelp6
            call add22(wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp7,fact)

          else if (sa3 == ssu) then

            ! def posB
            nhelp3 = b%d(ib,1)
            ! def fictive dimensions
            nhelp9 = nhelp6*(nhelp6-1)/2
            call add44(wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp9,nhelp6,fact)

          else
            ! sa3<ssu  B(pq,sr) <-- -A_s (pq,r)

            ! def posB-
            nhelp3 = b%d(ibm,1)
            call add32(wrk(nhelp2),wrk(nhelp3),u,nhelp8,nhelp7,nhelp6,-fact)

          end if

        end do

      else
        ! RC=16: nindB=4, nindext=1, typext=4, (typA is not 0 or 4, (NCI))
        rc = 16
        return
      end if

    else
      ! RC=17: nindB=4, nindext=1, typext=@  (Stup)
      rc = 17
      return
    end if

  else if (nindext == 2) then

    if (typext == 5) then

      if (typb == 0) then

        !4250 case B(p,q,r,s) <-- A(r,s)

        !     tests

        if ((typb == 0) .and. (typa /= 0)) then
          ! RC=18: nindB=4, nindext=2, typext=5, typB=0 (typA is not 0, Stup)
          rc = 18
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)

          nhelp1 = mmul(ssu,ssv)
          nhelp1 = mmul(nhelp1,sa1)
          nhelp1 = mmul(nhelp1,ssb)
          if (nhelp1 /= sa2) then
            write(u6,*) ' Add Bpqrs <- Ars incorrect',ssp,ssq,sa1,nhelp1,sa1,sa2
            cycle
          end if

          ib = b%i(ssu,ssv,sa1)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! calc joined pq index
          pq = (v-1)*nhelp4+u

          ! calc fictive lengths
          nhelp9 = nhelp6*nhelp7
          nhelp10 = nhelp4*nhelp5

          call add21(wrk(nhelp2),wrk(nhelp3),pq,nhelp10,nhelp9,fact)

        end do

      else if (typb == 4) then

        !4254 case B(pq,rs) <-- A(rs)

        !     tests

        if ((typb == 4) .and. (typa /= 1)) then
          ! RC=19: nindB=4, nindext=2, typext=5, typB=4 (typA is not 1, Stup)
          rc = 19
          return
        end if

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)

          nhelp1 = mmul(ssp,ssq)
          nhelp1 = mmul(nhelp1,sa1)
          nhelp1 = mmul(nhelp1,ssb)
          if (nhelp1 /= sa2) then
            write(u6,*) ' Add Bpqrs <- Ars incorrect'
            cycle
          end if

          ib = b%i(ssp,ssq,sa1)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr,dims
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))
          nhelp7 = dimm(b%d(0,4),b%d(ib,6))

          ! calc joined pq index and fictive length of pq pair
          if (ssp == ssq) then
            pq = (p-1)*(p-2)/2+q
            nhelp10 = nhelp4*(nhelp4-1)/2
          else
            pq = (q-1)*nhelp4+p
            nhelp10 = nhelp4*nhelp5
          end if

          ! calc fictive lengths
          if (sa1 == sa2) then
            nhelp9 = nhelp6*(nhelp6-1)/2
          else
            nhelp9 = nhelp6*nhelp7
          end if

          call add21(wrk(nhelp2),wrk(nhelp3),pq,nhelp10,nhelp9,fact)

        end do

      else
        ! RC=20: nindB=4, nindext=2, typext=5 (typB is not 0 or 4, NCI)
        rc = 20
        return
      end if

    else if (typext == 6) then

      !426  case B(p,q,r,s) <-- A(p,q) and B(pq,rs) <-- A(pq)

      ! RC=21: nindB=4, nindext=2, typext=6, NCI)
      rc = 21
      return

    else
      ! RC=22: nindB=4, nindext=2, (typext is not 5 or 6, NCI)
      rc = 22
      return
    end if

  else
    ! RC=23: nindB=4, nindext>2 (NCI)
    rc = 23
    return
  end if

else if (nindb == 3) then

  ! **********  -> B(pqr) **********

  if (nindext == 0) then

    !300  case B(pqr) <-- A(pqr)

    !     tests

    if (typa /= typb) then
      ! RC=24: nindB=3, nindext=0 (TypA incompatible with TypB ,Stup)
      rc = 24
      return
    end if

    do ia=1,a%d(0,5)

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)
      sa3 = a%d(ia,5)

      ib = b%i(sa1,sa2,1)

      ! def length
      nhelp1 = a%d(ia,2)
      if (nhelp1 == 0) cycle

      ! def posA,posB
      nhelp2 = a%d(ia,1)
      nhelp3 = b%d(ib,1)

      call add10(wrk(nhelp2),wrk(nhelp3),nhelp1,fact)

    end do

  else if (nindext == 1) then

    if (typext == 1) then

      !311  case B(p,q,r) <-- A(q,r)

      if ((typa == 0) .and. (typb == 0)) then

        !311  case B(p,q,r) <-- A(q,r)

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)

          ib = b%i(ssu,sa1,1)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))

          ! def fictive dimensions
          nhelp7 = nhelp5*nhelp6

          call add21(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp7,fact)

        end do

      else
        ! RC=25: nindB=3, nindext=1, typext=1 (tybA,B is not 0, NCI)
        rc = 25
        return
      end if

    else if (typext == 2) then

      if ((typa == 0) .and. (typb == 0)) then

        !312  case B(p,q,r) <-- A(p,r)

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)

          ib = b%i(sa1,ssu,1)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))

          call add32(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp6,nhelp7,fact)

        end do

      else
        ! RC=26: nindB=3, nindext=1, typext=2 (tybA,B is not 0, NCI)
        rc = 26
        return
      end if

    else if (typext == 3) then

      if ((typa == 0) .and. (typb == 0)) then

        !313  case B(p,q,r) <-- A(p,q)

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          sa2 = a%d(ia,4)

          ib = b%i(sa1,sa2,1)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq,dimr
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))
          nhelp6 = dimm(b%d(0,3),b%d(ib,5))

          ! def fictive dimensions
          nhelp7 = nhelp4*nhelp5

          call add21(wrk(nhelp2),wrk(nhelp3),u,nhelp7,nhelp6,fact)

        end do

      else
        ! RC=27: nindB=3, nindext=1, typext=3 (tybA,B is not 0, NCI)
        rc = 27
        return
      end if

    else
      ! RC=28: nindB=3 , typext=@ (Stup)
      rc = 28
      return
    end if

  else
    ! RC=29: nindB=3, nindext>1 (NCI)
    rc = 29
    return
  end if

else if (nindb == 2) then

  ! **********  -> B(pq) **********

  if (nindext == 0) then

    !200  case B(p,q) <-- A(p,q)

    do ia=1,a%d(0,5)

      sa1 = a%d(ia,3)
      sa2 = a%d(ia,4)

      ib = b%i(sa1,1,1)

      ! def length
      nhelp1 = a%d(ia,2)
      if (nhelp1 == 0) cycle

      ! def posA,posB
      nhelp2 = a%d(ia,1)
      nhelp3 = b%d(ib,1)

      call add10(wrk(nhelp2),wrk(nhelp3),nhelp1,fact)

    end do

  else if (nindext == 1) then

    if (typext == 1) then

      !211  case B(p,q) <-- A(q)

      if ((typa == 0) .and. (typb == 0)) then

        do ia=1,a%d(0,5)

          ib = b%i(ssu,1,1)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))

          call add21(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,fact)

        end do

      else
        ! RC=30: nindB=2, nindext=1, typext=1, (typA,B is not 0, NCI)
        rc = 30
        return
      end if

    else if (typext == 2) then

      !212  case B(p,q) <-- A(p)

      if ((typa == 0) .and. (typb == 0)) then

        do ia=1,a%d(0,5)

          sa1 = a%d(ia,3)
          ib = b%i(sa1,1,1)

          ! def length
          nhelp1 = a%d(ia,2)
          if (nhelp1 == 0) cycle

          ! def posA,posB
          nhelp2 = a%d(ia,1)
          nhelp3 = b%d(ib,1)

          ! def dimp,dimq
          nhelp4 = dimm(b%d(0,1),b%d(ib,3))
          nhelp5 = dimm(b%d(0,2),b%d(ib,4))

          call add22(wrk(nhelp2),wrk(nhelp3),u,nhelp4,nhelp5,fact)

        end do

      else
        ! RC=31: nindB=2, nindext=1, typext=2, (typA,B is not 0, NCI)
        rc = 31
        return
      end if

    else
      ! RC=32: nindB=2, nindext=1, typext=@ (Stup)
      rc = 32
      return
    end if

  else
    ! RC=33: nindB=2, ininext>1 (NCI)
    rc = 33
    return
  end if

else
  ! RC=34: nindb less then 2 (NCI/Stup)
  rc = 34
  return
end if

return

end subroutine add
