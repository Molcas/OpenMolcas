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

subroutine cct3_mult(wrk,wrksize,ninda,nindb,nindc,nindsum,a,ssa,b,ssb,c,ssc,rc)
! ninda   - # of indices in matrix A (Input)
! nindb   - # of indices in matrix B (Input)
! nindc   - # of indices in matrix C (Input, for test)
! nindsum - # of summation indices   (Input)
! a       - A  (Input)
! ssa     - overall symmetry state of matrix A (Input)
! b       - B  (Input)
! ssb     - overall symmetry state of matrix B (Input)
! c       - C  (Output)
! ssc     - overall symmetry state of matrix C (Input/Output)
!
! This routine realizes matrix-matrix and matrix-vector multiplications
! A(indA)*B(indB)=C(indC) or A(indA)*B(indB)=Y(indC)
!
! N.B.
! ninda, nindb, and nindsum are independent,
! nindc is tested if it is in agreement with previous ones - for test
! typA, typB are tested if they are in mutual agreement
!
! Table of implemented processes
!
! nindA   nindB  nindsum =>nindC        Operation               Implemented
!
! 4       4       4       0     A(pqrs)*B(pqrs)  =S             Not yet
! 4       4       3       2     A(p,qrs)*B(qrs,t)=C(p,t)          Yes
! 4       4       2       4     A(pq,rs)*B(rs,tu)=C(pq,tu)        Yes
! 4       4       1       6     A(pqr,s)*B(s,tuv)=C(pqr,tuv)      No
! 4       4       0       8     A(pqrs)*B(tuvz)=C(pqrs,tuvz)      No
!
! 4       3       3       1     A(p,qrs)*B(qrs)  =Y(p)            Yes
! 4       3       2       3     A(pq,rs)*B(rs,t) =C(pq,t)         Yes
! 4       3       1       5     A(pqr,s)*B(s,tu) =C(pqr,tu)       No
! 4       3       0       7     A(pqrs)*B(tuv)   =C(pqrs,tuv)     No
!
! 4       2       2       2     A(pq,rs)*B(pq)   =Y(pq)           Yes
! 4       2       1       4     A(pqr,s)*B(s,t)  =C(pqr,t)        Yes
! 4       2       0       6     A(pqrs)*B(tu)    =C(pqrs,tu)      No
!
! 4       1       1       3     A(pqr,s)*B(s)    =Y(pqr)        Not yet
! 4       1       0       5     A(pqrs)*B(t)     =C(pqrs,t)       No
!
! 3       4       3       1     A(pqr) *B(pqr,s) =Y(s)          Not need
! 3       4       2       3     A(p,qr)*B(qr,st) =C(p,st)         Yes
! 3       4       1       5     A(pq,r)*B(r,stu) =C(pq,stu)       No
! 3       4       0       7     A(pqr)*B(stuv) =C(pqr,stuv)       No
!
! 3       3       3       0     A(pqr) *B(pqr)   =S             Not yet
! 3       3       2       2     A(p,qr)*B(qr,s)  =C(p,s)          Yes
! 3       3       1       4     A(pq,r)*B(r,st)  =C(pq,st)        Yes
! 3       3       0       6     A(pqr) *B(stu)   =C(pqrmstu)      No
!
! 3       2       2       1     A(p,qr)*B(qr)    =Y(p)            Yes
! 3       2       1       3     A(pq,r)*B(r,s)   =C(pq,s)         Yes
! 3       2       0       5     A(pqr) *B(st)    =C(pqr,st)       No
!
! 3       1       1       2     A(pq,r)*B(r)     =Y(pq)         Not yet
! 3       1       0       4     A(pqr) *B(s)     =C(pqr,s)        No
!
! 2       4       2       2     A(pq) *B(pq,rs)  =Y(rs)         Not need
! 2       4       1       4     A(p,q)*B(q,rst)  =C(p,rst)        Yes
! 2       4       0       6     A(pq) *B(rstu)   =C(pq,rstu)      No
!
! 2       3       2       1     A(pq) *B(pq,r)   =Y(r)          Not need
! 2       3       1       3     A(p,q)*B(q,rs)   =C(p,rs)         Yes
! 2       3       0       5     A(pq) *B(rst)    =C(pq,rst)       No
!
! 2       2       2       0     A(pq) *B(pq)     =S             Not yet
! 2       2       1       2     A(p,q)*B(q,r)    =C(p,r)          Yes
! 2       2       0       4     A(pq) *B(rs)     =C(pq,rs)        No
!
! 2       1       1       1     A(p,q)*B(q)      =Y(p)          Not yet
! 2       1       0       3     A(pq) *B(r)      =C(pq,r)         No
!
! 1       4       1       3     A(p)*B(p,qrs)    =Y(qrs)        Not need
! 1       4       0       5     A(p)*B(qrst)     =C(p,qrst)       No
!
! 1       3       1       2     A(p)*B(p,qr)     =Y(qr)         Not need
! 1       3       0       4     A(p)*B(qrs)      =C(p,qrs)        No
!
! 1       2       1       1     A(p)*B(p,q)      =Y(q)          Not need
! 1       2       0       3     A(p)*B(qr)       =C(p,qr)         No
!
! 1       1       1       0     A(p)*B(p)        =S             Not yet
! 1       1       0       2     A(p)*B(q)        =C(p,q)          No
!
! Legend used in error messages:
! NCI = Not Currently Implemented
! Stup= Stupidity
! @   = improper value

use CCT3_global, only: Map_Type, mmul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, ninda, nindb, nindc, nindsum, ssa, ssb
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: ssc, rc
integer(kind=iwp) :: ix, mvec(4096,7), typa, typb

rc = 0
ssc = mmul(ssa,ssb)
typa = a%d(0,6)
typb = b%d(0,6)

if (ninda == 4) then

  ! ****** nind A = 4  ******

  if (nindb == 4) then

    !44 *** nind A =4, B=4

    if (nindsum == 1) then
      ! RC=1 : nindA=4, nindB=4, nindsum=1 (too large metiate, NCI)
      rc = 1
      return

    else if (nindsum == 2) then

      ! case A(pqrs)*B(rstu)=C(pqtu)

      ! tests

      if (nindc /= 4) then
        ! RC=2 : nindA=4, nindB=4, nindsum=2 (# of ind. in C is not 4, Stup)
        rc = 2
        return
      end if

      if (typa == 2) then
        ! RC=3 : nindA=4, nindB=4, nindsum=2 (typA is 2, Stup)
        rc = 3
        return
      end if

      if (typb == 2) then
        ! RC=4 : nindA=4, nindB=4, nindsum=2 (typB is 2, Stup)
        rc = 4
        return
      end if

      if ((typa == 3) .or. (typa == 4)) then
        ! there is  r>s in A
        if ((typb == 4) .or. (typb == 1)) then
          ! OK
        else
          ! RC=5 : nindA=4, nindB=4, nindsum=2 (typA incompatible with typB , Stup)
          rc = 5
          return
        end if
      else
        ! there is no r>s in A
        if ((typb == 4) .or. (typb == 1)) then
          ! RC=5 : nindA=4, nindB=4, nindsum=2 (typA incompatible with typB , Stup)
          rc = 5
          return
        else
          ! OK
        end if
      end if

      ! call cct3_grc44c and multc0

      call cct3_grc44c(a,b,c,mvec,ssa,ssb,2,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 3) then

      ! case A(pqrs)*B(qrst)=C(pt)

      ! tests

      if (nindc /= 2) then
        ! RC=6 : nindA=4, nindB=4, nindsum=3 (# of ind. in C is not 2, Stup)
        rc = 6
        return
      end if

      if (typa == 1) then
        ! RC=7 : nindA=4, nindB=4, nindsum=3 (typA is 1, Stup)
        rc = 7
        return
      end if

      if (typb == 3) then
        ! RC=8 : nindA=4, nindB=4, nindsum=3 (typB is 3, Stup)
        rc = 8
        return
      end if

      if ((typa == 2) .and. (typb /= 1)) then
        ! RC=9 : nindA=4, nindB=4, nindsum=3 (typA incomp. with typB, Stup)
        rc = 9
        return
      end if

      if ((typa == 3) .and. (typb /= 2)) then
        ! RC=9 : nindA=4, nindB=4, nindsum=3 (typA incomp. with typB, Stup)
        ! like in previous case
        rc = 9
        return
      end if

      ! call cct3_grc44c and multc0

      call cct3_grc44c(a,b,c,mvec,ssa,ssb,1,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 4) then
      ! RC=10: nindA=4, nindB=4, nindsum=4 (NCI)
      rc = 10
      return

    else
      ! RC=11: nindA=4, nindB=4, nindsum=@ (Stup)
      rc = 11
      return
    end if

  else if (nindb == 3) then

    !43 *** nind A =4, B=3

    if (nindsum == 3) then

      ! case A(pqrs)*B(qrs)=Y(p)

      ! tests

      if (nindc /= 1) then
        ! RC=12: nindA=4, nindB=3, nindsum=3 (# of index in C is not 1, Stup)
        rc = 12
        return
      end if

      if ((typa == 1) .or. (typa == 4)) then
        ! RC=13: nindA=4, nindB=3, nindsum=3 (typA is 1 or 4, Stup)
        rc = 13
        return
      end if

      if ((typa == 2) .and. (typb /= 1)) then
        ! RC=14: nindA=4, nindB=3, nindsum=3 (typA incomp. with typB, Stup)
        rc = 14
        return
      end if

      if ((typa == 3) .and. (typb /= 2)) then
        ! RC=14: nindA=4, nindB=3, nindsum=3 (typA incomp. with typB, Stup)
        ! like in previous case
        rc = 14
        return
      end if

      ! call cct3_grc43y and multy0

      call cct3_grc43y(a,b,c,mvec,ssa,ssb,ix)

      call cct3_multy0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 2) then

      ! case A(pqrs)*B(rst)=C(pqt)

      ! tests

      if (nindc /= 3) then
        ! RC=15: nindA=4, nindB=3, nindsum=2 (# of index in C is not 3, Stup)
        rc = 15
        return
      end if

      if (typa == 2) then
        ! RC=16: nindA=4, nindB=3, nindsum=2 (typA is 2, Stup)
        rc = 16
        return
      end if

      if (typb == 2) then
        ! RC=17: nindA=4, nindB=3, nindsum=2 (typB is 2, Stup)
        rc = 17
        return
      end if

      if (((typa == 3) .or. (typa == 4)) .and. (typb /= 1)) then
        ! RC=18: nindA=4, nindB=3, nindsum=2 (typA incomp. with typB, Stup)
        rc = 18
        return
      end if

      ! call cct3_grc43c and multc0

      call cct3_grc43c(a,b,c,mvec,ssa,ssb,2,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 1) then
      ! RC=19: nindA=4, nindB=3, nindsum=1 (too large mediate, NCI)
      rc = 19
      return

    else
      ! RC=20: nindA=4, nindB=3, nindsum=@ (Stup)
      rc = 20
      return
    end if

  else if (nindb == 2) then

    !42 *** nind A =4, B=2

    if (nindsum == 1) then

      ! case A(pqr,s)*B(s,t)=C(pqr,t)

      ! tests

      if (nindc /= 4) then
        ! RC=21: nindA=4, nindB=2, nindsum=1 (# of index in C is not 4, Stup)
        rc = 21
        return
      end if

      if ((typa == 3) .or. (typa == 4)) then
        ! RC=22: nindA=4, nindB=2, nindsum=1 (typA is 3 or 4, Stup)
        rc = 22
        return
      end if

      if (typb == 1) then
        ! RC=23: nindA=4, nindB=2, nindsum=1 (typA is 1, Stup)
        rc = 23
        return
      end if

      ! call cct3_grc42c and multc0

      call cct3_grc42c(a,b,c,mvec,ssa,ssb,3,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 2) then

      ! case A(pqrs)*B(rs)=Y(pq)

      ! tests

      if (nindc /= 2) then
        ! RC=24: nindA=4, nindB=2, nindsum=2 (# of index in C is not 2, Stup)
        rc = 24
        return
      end if

      if (typa == 2) then
        ! RC=25: nindA=4, nindB=2, nindsum=2 (typA is 2, Stup)
        rc = 25
        return
      end if

      if (((typa == 3) .or. (typa == 4)) .and. (typb /= 1)) then
        ! RC=26: nindA=4, nindB=2, nindsum=2 (typA incomp. with typB, Stup)
        rc = 26
        return
      end if

      ! call cct3_grc42y and multy0

      call cct3_grc42y(a,b,c,mvec,ssa,ssb,ix)

      call cct3_multy0(wrk,wrksize,mvec,ix,c,1)

    else
      ! RC=27: nindA=4, nindB=2, nindsum=@ (Stup)
      rc = 27
      return
    end if

  else if (nindb == 1) then

    !41 *** nind A =4, B=1

    ! RC=28: nindA=4, nindB=1 (NCI)
    rc = 28
    return

  else

    !4@ *** nind A =4, B=@

    ! RC=29  : nindA=4, nindB=@ (Stup)
    rc = 29
    return
  end if

else if (ninda == 3) then

  ! ****** nind A = 3  ******

  if (nindb == 4) then

    !34 *** nind A =3, B=4

    if (nindsum == 1) then
      ! RC=30: nindA=3, nindB=4, nindsum=1 (too large mediate, NCI)
      rc = 30
      return

    else if (nindsum == 2) then

      ! case A(p,qr)*B(qr,st)=C(p,st)

      ! tests

      if (nindc /= 3) then
        ! RC=31: nindA=3, nindB=4, nindsum=2 (# of ind. in C is not 3, Stup)
        rc = 31
        return
      end if

      if (typa == 1) then
        ! RC=32: nindA=3, nindB=4, nindsum=2 (typA is 1, Stup)
        rc = 32
        return
      end if

      if (typb == 2) then
        ! RC=33: nindA=3, nindB=4, nindsum=2 (typB is 2, Stup)
        rc = 33
        return
      end if

      if (((typb == 1) .or. (typb == 4)) .and. (typa /= 2)) then
        ! RC=34: nindA=3, nindB=4, nindsum=2 (typA incomp. with typB, Stup)
        rc = 34
        return
      end if

      ! call cct3_grc34c and multc0

      call cct3_grc34c(a,b,c,mvec,ssa,ssb,1,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 3) then
      ! RC=35: nindA=3, nindB=4, nindsum=3 (NCI)
      rc = 35
      return
    end if

  else if (nindb == 3) then

    !33 *** nind A =3, B=3

    if (nindsum == 1) then

      ! case A(pq,r)*B(r,st)=A(pq,st)

      ! tests

      if (nindc /= 4) then
        ! RC=36: nindA=3, nindB=3, nindsum=1 (# of ind. in C is not 4, Stup)
        rc = 36
        return
      end if

      if (typa == 2) then
        ! RC=37: nindA=3, nindB=3, nindsum=1 (typA is 2, Stup)
        rc = 37
        return
      end if

      if (typb == 1) then
        ! RC=38: nindA=3, nindB=3, nindsum=1 (typB is 1, Stup)
        rc = 38
        return
      end if

      ! call cct3_grc34c and multc0

      call cct3_grc34c(a,b,c,mvec,ssa,ssb,2,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 2) then

      ! case A(p,qr)*B(rq,s)=C(p,s)

      ! tests

      if (nindc /= 2) then
        ! RC=39: nindA=3, nindB=3, nindsum=2 (# of ind. in C is not 2, Stup)
        rc = 39
        return
      end if

      if (typa == 1) then
        ! RC=40: nindA=3, nindB=3, nindsum=2 (typA is 1, Stup)
        rc = 40
        return
      end if

      if (typb == 2) then
        ! RC=41: nindA=3, nindB=3, nindsum=2 (typB is 2, Stup)
        rc = 41
        return
      end if

      if ((typa == 2) .and. (typb /= 1)) then
        ! RC=42: nindA=3, nindB=3, nindsum=2 (typA incomp. with typB, Stup)
        rc = 42
        return
      end if

      ! call cct3_grc33c and multc0

      call cct3_grc33c(a,b,c,mvec,ssa,ssb,1,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 3) then
      ! RC=43: nindA=3, nindB=3, nindsum=3 (NCI)
      rc = 43
      return
    end if

  else if (nindb == 2) then

    !32 *** nind A =3, B=2

    if (nindsum == 1) then

      ! case A(pq,r)*B(r,s)=C(pq,s)

      ! tests

      if (nindc /= 3) then
        ! RC=44: nindA=3, nindB=2, nindsum=1 (# of ind. in C is not 3, Stup)
        rc = 44
        return
      end if

      if (typa == 2) then
        ! RC=45: nindA=3, nindB=2, nindsum=1 (typA is 2, Stup)
        rc = 45
        return
      end if

      if (typb == 1) then
        ! RC=46: nindA=3, nindB=2, nindsum=1 (typB is 1, Stup)
        rc = 46
        return
      end if

      ! call cct3_grc32c and multc0

      call cct3_grc32c(a,b,c,mvec,ssa,ssb,2,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 2) then

      ! case A(p,qr)*B(qr)=C(p)

      ! tests

      if (nindc /= 1) then
        ! RC=47: nindA=3, nindB=2, nindsum=2 (# of ind. in C is not 1, Stup)
        rc = 47
        return
      end if

      if (typa == 1) then
        ! RC=48: nindA=3, nindB=2, nindsum=2 (typA is 1, Stup)
        rc = 48
        return
      end if

      if ((typb == 1) .and. (typa /= 2)) then
        ! RC=49: nindA=3, nindB=2, nindsum=2 (typA incomp. with typB, Stup)
        rc = 49
        return
      end if

      ! call cct3_grc32y and multy0

      call cct3_grc32y(a,b,c,mvec,ssa,ssb,ix)

      call cct3_multy0(wrk,wrksize,mvec,ix,c,1)

    else
      ! RC=50: nindA=3, nindB=2, nindsum=@ (Stup)
      rc = 50
      return
    end if

  else if (nindb == 1) then

    !31 *** nind A =3, B=1

    ! RC=51: nindA=3, nindB=1 (NCI)
    rc = 51
    return

  else

    !3@ *** nind A =3, B=@

    ! RC=52: nindA=3, nindB=@ (Stup)
    rc = 52
    return
  end if

else if (ninda == 2) then

  ! ****** nind A = 2  ******

  if (nindb == 4) then

    !24 *** nind A =2, B=4

    if (nindsum == 1) then

      ! case A(p,q)*B(q,rst)=C(p,rst)

      ! tests

      if (nindc /= 4) then
        ! RC=53: nindA=2, nindB=4, nindsum=1 (# of ind. in C is not 4, Stup)
        rc = 53
        return
      end if

      if (typa == 1) then
        ! RC=54: nindA=2, nindB=4, nindsum=1 (typA is 1, Stup)
        rc = 54
        return
      end if

      if ((typb == 1) .or. (typb == 4)) then
        ! RC=55: nindA=2, nindB=4, nindsum=1 (typB is 1 or 4, Stup)
        rc = 55
        return
      end if

      ! call cct3_grc24c and multc0

      call cct3_grc24c(a,b,c,mvec,ssa,ssb,1,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 2) then
      ! RC=56: nindA=2, nindB=4, nindsum=2 (NCI)
      rc = 56
      return

    else
      ! RC=57: nindA=2, nindB=4, nindsum=@ (Stup)
      rc = 57
      return
    end if

  else if (nindb == 3) then

    !23 *** nind A =2, B=3

    if (nindsum == 1) then

      ! case A(p,q)*B(q,rs)=C(p,rs)

      ! tests

      if (nindc /= 3) then
        ! RC=58: nindA=2, nindB=3, nindsum=1 (# of ind. in C is not 3, Stup)
        rc = 58
        return
      end if

      if (typa == 1) then
        ! RC=59: nindA=2, nindB=3, nindsum=1 (typA is 1, Stup)
        rc = 59
        return
      end if

      if (typb == 1) then
        ! RC=60: nindA=2, nindB=3, nindsum=1 (typB is 1, Stup)
        rc = 60
        return
      end if

      ! call cct3_grc23c and multc0

      call cct3_grc23c(a,b,c,mvec,ssa,ssb,1,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 2) then
      ! RC=61: nindA=2, nindB=3, nindsum=2 (NCI)
      rc = 61
      return

    else
      ! RC=62: nindA=2, nindB=3, nindsum=@ (Stup)
      rc = 62
      return
    end if

  else if (nindb == 2) then

    !22 *** nind A =2, B=2

    if (nindsum == 1) then

      ! case A(p,q)*B(q,r)=C(p,r)

      ! tests

      if (nindc /= 2) then
        ! RC=63: nindA=2, nindB=2, nindsum=1 (# of ind. in C is not 2, Stup)
        rc = 63
        return
      end if

      if (typa == 1) then
        ! RC=64: nindA=2, nindB=2, nindsum=1 (typA is 1, Stup)
        rc = 64
        return
      end if

      if (typb == 1) then
        ! RC=65: nindA=2, nindB=2, nindsum=1 (typB is 1, Stup)
        rc = 65
        return
      end if

      ! call cct3_grc22c and multc0

      call cct3_grc22c(a,b,c,mvec,ssa,ssb,1,ix)

      call cct3_multc0(wrk,wrksize,mvec,ix,c,1)

    else if (nindsum == 2) then
      ! RC=66: nindA=2, nindB=2, nindsum=2 (NCI)
      rc = 66
      return

    else
      ! RC=67: nindA=2, nindB=2, nindsum=@ (Stup)
      rc = 67
      return
    end if

  else if (nindb == 1) then

    !21 *** nind A =2, B=1

    ! RC=68: nindA=2, nindB=1 (NCI)
    rc = 68
    return

  else

    !2@ *** nind A =2, B=@

    ! RC=69: nindA=2, nindB=@
    rc = 69
    return
  end if

else if (ninda == 1) then

  ! ****** nind A = 1  ******

  !1x RC=70: nind=1 (NCI)
  rc = 70
  return

else

  ! ****** nind A = @  ******

  !@x RC=71: inproper nindA
  rc = 71
  return
end if

return

end subroutine cct3_mult
