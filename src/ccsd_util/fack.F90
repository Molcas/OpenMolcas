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

subroutine fack(wrk,wrksize,nind,newtyp,mapda,ssa,mapia,mapdb,mapib,possb0,rc)
! nind   - # of indices in matrices A,B (Input)
! newtyp - typ of final matrix B (Input) (i.e. mapdb(0,6))
! mapda  - direct map matrix corresponding to A (Input)
! ssa    - overall symmetry state of matrix A (Input)
! mapia  - inverse map matrix corresponding to A (Input)
! mapdb  - direct map matrix corresponding to B (Output)
! mapib  - inverse map matrix corresponding to B (Output)
! possb0 - initial position of B matrix in WRK (Input)
! rc     - return (error) code (Output)
!
! this routine realizes following packings :
!
! Operation                               Nind    TypA    TypB
!
! B(pq)     = A(p,q)     - A(q,p)          2       0       1
! B(pq,r)   = A(p,q,r)   - A(q,p,r)        3       0       1
! B(p,qr)   = A(p,q,r)   - A(p,r,q)        3       0       2
! B(pq,r,s) = A(p,q,r,s) - A(q,p,r,s)      4       0       1
! B(p,q,rs) = A(p,q,r,s) - A(q,p,s,r)      4       0       3
! B(pq,rs)  = NCI                          4       0       4
! B(pq,rs)  = A(pq,r,s)  - A(pq,s,r)       4       1       4
! B(pq,rs)  = A(p,q,rs)  - A(q,p,rs)       4       3       4

#include "ccsd1.fh"
#include "wrk.fh"
integer nind, newtyp, ssa, possb0, rc
integer mapda(0:512,1:6)
integer mapdb(0:512,1:6)
integer mapia(1:8,1:8,1:8)
integer mapib(1:8,1:8,1:8)
! help variables
integer nhelp1, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6, nhelp7, nhelp8, nhelp9
integer sb1, sb2, sb3, sb4
integer rc1, ib, iap, iam, posst

rc = 0

! get mapdb,mapib

call grc0(nind,newtyp,mapda(0,1),mapda(0,2),mapda(0,3),mapda(0,4),ssa,possb0,posst,mapdb,mapib)

if (nind < 2) then

  ! *********  less than 2 ind ***

  ! RC=1 less than 2 indices
  rc = 1
  return

else if (nind == 2) then

  ! *********  2 indices *********

  !2.1 A(p,q) -> B(pq)

  if (mapda(0,6) /= 0) then
    ! RC=2 bad type for: nind=2
    rc = 2
    return
  end if

  do ib=1,mapdb(0,5)

    sb1 = mapdb(ib,3)
    sb2 = mapdb(ib,4)
    if (mapdb(ib,2) == 0) cycle

    if (sb1 == sb2) then
      ! sym b1 = sym b2

      ! def ia+,-
      iap = mapia(sb1,1,1)

      ! poss B,A+,A-
      nhelp1 = mapdb(ib,1)
      nhelp2 = mapda(iap,1)

      ! dimp-s
      nhelp4 = dimm(mapdb(0,1),sb1)

      ! def size
      nhelp9 = nhelp4*(nhelp4-1)/2

      call pack210(wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp4,rc1)

    else
      ! sym b1 > sym b2

      ! def ia+,-
      iap = mapia(sb1,1,1)
      iam = mapia(sb2,1,1)

      ! poss B,A+,A-
      nhelp1 = mapdb(ib,1)
      nhelp2 = mapda(iap,1)
      nhelp3 = mapda(iam,1)

      ! dimp-s
      nhelp4 = dimm(mapdb(0,1),sb1)
      nhelp5 = dimm(mapdb(0,2),sb2)

      call pack211(wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,rc1)

    end if

  end do

else if (nind == 3) then

  ! *********  3 indices *********

  if (mapda(0,6) /= 0) then
    ! RC=3 bad type for: nind=3
    rc = 3
    return
  end if

  if (newtyp == 1) then

    !3.1 A(p,q,r) -> B(pq,r)

    do ib=1,mapdb(0,5)

      sb1 = mapdb(ib,3)
      sb2 = mapdb(ib,4)
      sb3 = mapdb(ib,5)
      if (mapdb(ib,2) == 0) cycle

      if (sb1 == sb2) then
        ! sym b1 = sym b2

        ! def ia+,-
        iap = mapia(sb1,sb2,1)

        ! poss B,A+,A-
        nhelp1 = mapdb(ib,1)
        nhelp2 = mapda(iap,1)

        ! dimp-s
        nhelp4 = dimm(mapdb(0,1),sb1)
        nhelp6 = dimm(mapdb(0,3),sb3)

        ! def size
        nhelp9 = nhelp4*(nhelp4-1)/2

        call pack310(wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp6,nhelp4,rc1)

      else
        ! sym b1 > sym b2

        ! def ia+,-
        iap = mapia(sb1,sb2,1)
        iam = mapia(sb2,sb1,1)

        ! poss B,A+,A-
        nhelp1 = mapdb(ib,1)
        nhelp2 = mapda(iap,1)
        nhelp3 = mapda(iam,1)

        ! dimp-s
        nhelp4 = dimm(mapdb(0,1),sb1)
        nhelp5 = dimm(mapdb(0,2),sb2)
        nhelp6 = dimm(mapdb(0,3),sb3)

        call pack311(wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,nhelp6,rc1)

      end if

    end do

  else if (newtyp == 2) then

    !3.2 A(p,q,r) -> B(p,qr)

    do ib=1,mapdb(0,5)

      sb1 = mapdb(ib,3)
      sb2 = mapdb(ib,4)
      sb3 = mapdb(ib,5)
      if (mapdb(ib,2) == 0) cycle

      if (sb2 == sb3) then
        ! sym b2 = sym b3

        ! def ia+,-
        iap = mapia(sb1,sb2,1)

        ! poss B,A+,A-
        nhelp1 = mapdb(ib,1)
        nhelp2 = mapda(iap,1)

        ! dimp-s
        nhelp4 = dimm(mapdb(0,1),sb1)
        nhelp5 = dimm(mapdb(0,2),sb2)

        ! def size
        nhelp9 = nhelp5*(nhelp5-1)/2

        call pack320(wrk(nhelp2),wrk(nhelp1),nhelp4,nhelp9,nhelp5,rc1)

      else
        ! sym b2 > sym b3

        ! def ia+,-
        iap = mapia(sb1,sb2,1)
        iam = mapia(sb1,sb3,1)

        ! poss B,A+,A-
        nhelp1 = mapdb(ib,1)
        nhelp2 = mapda(iap,1)
        nhelp3 = mapda(iam,1)

        ! dimp-s
        nhelp4 = dimm(mapdb(0,1),sb1)
        nhelp5 = dimm(mapdb(0,2),sb2)
        nhelp6 = dimm(mapdb(0,3),sb3)

        call pack321(wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,nhelp6,rc1)

      end if

    end do

  else
    ! RC=4 : bad newtyp for: nind=3
    rc = 4
    return
  end if

else if (nind == 4) then

  ! *********  4 indices *********

  if (mapda(0,6) == 0) then

    !* case A(p,q,r,s) => newtyp can be 1,3,4

    if (newtyp == 1) then

      !4.1 A(p,q,r,s) -> B(pq,r,s)

      do ib=1,mapdb(0,5)

        sb1 = mapdb(ib,3)
        sb2 = mapdb(ib,4)
        sb3 = mapdb(ib,5)
        sb4 = mapdb(ib,6)
        if (mapdb(ib,2) == 0) cycle

        if (sb1 == sb2) then
          ! sym b1 = sym b2

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp6 = dimm(mapdb(0,3),sb3)
          nhelp7 = dimm(mapdb(0,4),sb4)

          ! def size
          nhelp8 = nhelp4*(nhelp4-1)/2
          nhelp9 = nhelp6*nhelp7

          call pack310(wrk(nhelp2),wrk(nhelp1),nhelp8,nhelp9,nhelp4,rc1)

        else
          ! sym b1 > sym b2

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)
          iam = mapia(sb2,sb1,sb3)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)
          nhelp3 = mapda(iam,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp5 = dimm(mapdb(0,2),sb2)
          nhelp6 = dimm(mapdb(0,3),sb3)
          nhelp7 = dimm(mapdb(0,4),sb4)

          ! def size
          nhelp9 = nhelp6*nhelp7

          call pack311(wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,nhelp9,rc1)

        end if

      end do

    else if (newtyp == 3) then

      !4.2 A(p,q,r,s) -> B(p,q,rs)

      do ib=1,mapdb(0,5)

        sb1 = mapdb(ib,3)
        sb2 = mapdb(ib,4)
        sb3 = mapdb(ib,5)
        sb4 = mapdb(ib,6)
        if (mapdb(ib,2) == 0) cycle

        if (sb3 == sb4) then
          ! sym b3 = sym b4

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp5 = dimm(mapdb(0,2),sb2)
          nhelp6 = dimm(mapdb(0,3),sb3)

          ! def size
          nhelp8 = nhelp6*(nhelp6-1)/2
          nhelp9 = nhelp4*nhelp5

          call pack320(wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp8,nhelp6,rc1)

        else
          ! sym b3 > sym b4

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)
          iam = mapia(sb1,sb2,sb4)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)
          nhelp3 = mapda(iam,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp5 = dimm(mapdb(0,2),sb2)
          nhelp6 = dimm(mapdb(0,3),sb3)
          nhelp7 = dimm(mapdb(0,4),sb4)

          ! def size
          nhelp9 = nhelp4*nhelp5

          call pack321(wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp9,nhelp6,nhelp7,rc1)

        end if

      end do

    else if (newtyp == 4) then

      !4.3 A(p,q,r,s) -> B(pq,rs)

      ! RC=5 : Not Currently Implemented : nind=4, oldtyp=0, newtyp=4
      rc = 5
      return

    else
      ! RC=6 : incompatible newtyp for: nind=4, oldtyp=0
      rc = 6
      return
    end if

  else if (mapda(0,6) == 1) then

    !* case A(pq,r,s) => newtyp can be 4

    if (newtyp == 4) then

      !4.4 A(pq,r,s) -> B(pq,rs)

      do ib=1,mapdb(0,5)

        sb1 = mapdb(ib,3)
        sb2 = mapdb(ib,4)
        sb3 = mapdb(ib,5)
        sb4 = mapdb(ib,6)
        if (mapdb(ib,2) == 0) cycle

        if (sb3 == sb4) then
          ! sym b3 = sym b4

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp5 = dimm(mapdb(0,2),sb2)
          nhelp6 = dimm(mapdb(0,3),sb3)

          ! def size
          nhelp8 = nhelp6*(nhelp6-1)/2
          if (sb1 == sb2) then
            nhelp9 = nhelp4*(nhelp4-1)/2
          else
            nhelp9 = nhelp4*nhelp5
          end if

          call pack320(wrk(nhelp2),wrk(nhelp1),nhelp9,nhelp8,nhelp6,rc1)

        else
          ! sym b3 > sym b4

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)
          iam = mapia(sb1,sb2,sb4)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)
          nhelp3 = mapda(iam,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp5 = dimm(mapdb(0,2),sb2)
          nhelp6 = dimm(mapdb(0,3),sb3)
          nhelp7 = dimm(mapdb(0,4),sb4)

          ! def size
          if (sb1 == sb2) then
            nhelp9 = nhelp4*(nhelp4-1)/2
          else
            nhelp9 = nhelp4*nhelp5
          end if

          call pack321(wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp9,nhelp6,nhelp7,rc1)

        end if

      end do

    else
      ! RC=7 : incompatible newtyp for: nind=4, oldtyp=1
      rc = 7
      return
    end if

  else if (mapda(0,6) == 3) then

    !* case A(p,q,rs) => newtyp can be 4

    if (newtyp == 4) then

      !4.5 A(p,q,rs) -> B(pq,rs)

      do ib=1,mapdb(0,5)

        sb1 = mapdb(ib,3)
        sb2 = mapdb(ib,4)
        sb3 = mapdb(ib,5)
        sb4 = mapdb(ib,6)
        if (mapdb(ib,2) == 0) cycle

        if (sb1 == sb2) then
          ! sym b1 = sym b2

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp6 = dimm(mapdb(0,3),sb3)
          nhelp7 = dimm(mapdb(0,4),sb4)

          ! def size
          nhelp8 = nhelp4*(nhelp4-1)/2
          if (sb3 == sb4) then
            nhelp9 = nhelp6*(nhelp6-1)/2
          else
            nhelp9 = nhelp6*nhelp7
          end if

          call pack310(wrk(nhelp2),wrk(nhelp1),nhelp8,nhelp9,nhelp4,rc1)

        else
          ! sym b1 > sym b2

          ! def ia+,-
          iap = mapia(sb1,sb2,sb3)
          iam = mapia(sb2,sb1,sb3)

          ! poss B,A+,A-
          nhelp1 = mapdb(ib,1)
          nhelp2 = mapda(iap,1)
          nhelp3 = mapda(iam,1)

          ! dimp-s
          nhelp4 = dimm(mapdb(0,1),sb1)
          nhelp5 = dimm(mapdb(0,2),sb2)
          nhelp6 = dimm(mapdb(0,3),sb3)
          nhelp7 = dimm(mapdb(0,4),sb4)

          ! def size
          if (sb3 == sb4) then
            nhelp9 = nhelp6*(nhelp6-1)/2
          else
            nhelp9 = nhelp6*nhelp7
          end if

          call pack311(wrk(nhelp2),wrk(nhelp3),wrk(nhelp1),nhelp4,nhelp5,nhelp9,rc1)

        end if

      end do

    else
      ! RC=8 : incompatible newtyp for: nind=4, oldtyp=3
      rc = 8
      return
    end if

  else
    ! RC=9 : incompatible typold for: nind=4
    rc = 9
    return
  end if

else

  ! *********  more than 4 ind ***

  ! RC=10: more than 4 indices
  rc = 10
  return
end if

return

end subroutine fack
