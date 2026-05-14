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

subroutine ccsort_t3grc0(nind,typ,typp,typq,typr,typs,stot,post,map)
! N.B. This routine is in principle copy of those in T3,
! but some changes were done:
! 1) mmul is substituted by mul
! 2) dimm is added, since using of ccsd1.com is impossible
! 3) ccsd1.com is replaced by ccsort_global
!
! nind   - number of indices (I)
! typ    - typ of mediate (I)
! typp   - typ of index p (I)
! typq   - typ of index q (I)
! typr   - typ of index r (I)
! typs   - typ of index s (I)
! stot   - overall symmetry of the mediate (I)
! post   - final position of the mediate (O)
! map    - map type of the mediate (O)
!
! this routine defines %d and %i for given intermediate
! it can done exactly the same maps like grc0 in CCSD
! plus additional types of mediates are introduced:
! type    meaning
! 5       p>q>r,s ; also p>q>r
! 6       p,q>r>s
! 7       p>=q,r,s ; also p>=q,r; p>=q
! 8       p,q>=r,s ; also p,q>=s
! 9       p,q,q>=s
! 10      p>=q,r>=s
! 11      p>=q>=r,s ; also p>=q>=r
! 12      p,q>=r>=s
!
! currently, these new types are implemented only for nind=3
!
! N.B. (this routine cannot run with +OP2)
! N.B. this routine does not test stupidities

use ccsort_global, only: Map_Type, noa, nob, NSYM, nva, nvb
use Symmetry_Info, only: Mul
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nind, typ, typp, typq, typr, typs, stot
integer(kind=iwp), intent(out) :: post
type(Map_Type), intent(_OUT_) :: map
integer(kind=iwp) :: dimm(5,8), i, nhelp1, nhelp2, nhelp3, nhelp4, nsymq, nsymr, pos, rsk1, rsk2, sp, spq, spqr, sq, sr, ss

! !!!!!!!! def dimm to je tu len terazky, lebo nemozeme pouzivat ccsd1.com !!!!

! Tutok musim cosi inicializovat
ss = 0
pos = 0
rsk1 = 0
rsk2 = 0
dimm(1,1:nsym) = noa(1:nsym)
dimm(2,1:nsym) = nob(1:nsym)
dimm(3,1:nsym) = nva(1:nsym)
dimm(4,1:nsym) = nvb(1:nsym)
dimm(5,1:nsym) = nva(1:nsym)+noa(1:nsym)

! vanishing %i files

map%i(1:nsym,1:nsym,1:nsym) = 0

i = 1
pos = map%pos0
if (nind == 1) then

  ! matrix A(p)

  sp = mul(stot,1)

  nhelp1 = dimm(typp,sp)

  ! def map%i
  map%i(1,1,1) = i

  ! def position
  map%d(i,1) = pos

  ! def length
  map%d(i,2) = nhelp1

  ! def sym p,q
  map%d(i,3) = sp
  map%d(i,4) = 0
  map%d(i,5) = 0
  map%d(i,6) = 0

  pos = pos+map%d(i,2)
  i = i+1

else if (nind == 2) then

  ! matrix A(p,q)

  do sp=1,nsym

    sq = mul(stot,sp)
    ! Meggie out
    if ((typ == 1) .and. (sp < sq)) cycle

    nhelp1 = dimm(typp,sp)
    nhelp2 = dimm(typq,sq)

    ! def map%i
    map%i(sp,1,1) = i

    ! def position
    map%d(i,1) = pos

    ! def length
    if ((typ == 1) .and. (sp == sq)) then
      map%d(i,2) = nhelp1*(nhelp1-1)/2
    else
      map%d(i,2) = nhelp1*nhelp2
    end if

    ! def sym p,q
    map%d(i,3) = sp
    map%d(i,4) = sq
    map%d(i,5) = 0
    map%d(i,6) = 0

    pos = pos+map%d(i,2)
    i = i+1

  end do

else if (nind == 3) then

  ! matrix A(p,q,r)

  ! def reucion summation keys : rsk1 for pq, rsk2 for qr

  if (typ == 0) then
    rsk1 = 0
    rsk2 = 0
  else if (typ == 1) then
    rsk1 = 1
    rsk2 = 0
  else if (typ == 2) then
    rsk1 = 0
    rsk2 = 1
  else if (typ == 5) then
    rsk1 = 1
    rsk2 = 1
  else if (typ == 7) then
    rsk1 = 1
    rsk2 = 0
  else if (typ == 8) then
    rsk1 = 0
    rsk2 = 1
  else if (typ == 11) then
    rsk1 = 1
    rsk2 = 1
  end if

  do sp=1,nsym
    if (rsk1 == 1) then
      nsymq = sp
    else
      nsymq = nsym
    end if

    do sq=1,nsymq
      spq = mul(sp,sq)

      sr = mul(stot,spq)
      ! Meggie out
      if ((rsk2 == 1) .and. (sq < sr)) cycle

      nhelp1 = dimm(typp,sp)
      nhelp2 = dimm(typq,sq)
      nhelp3 = dimm(typr,sr)

      ! def map%i
      map%i(sp,sq,1) = i

      ! def position
      map%d(i,1) = pos

      ! def length
      if ((typ == 1) .and. (sp == sq)) then
        map%d(i,2) = nhelp1*(nhelp1-1)*nhelp3/2
      else if ((typ == 2) .and. (sq == sr)) then
        map%d(i,2) = nhelp1*nhelp2*(nhelp2-1)/2
      else if (typ == 5) then
        if (sp == sr) then
          map%d(i,2) = nhelp1*(nhelp1-1)*(nhelp1-2)/6
        else if (sp == sq) then
          map%d(i,2) = nhelp1*(nhelp1-1)*nhelp3/2
        else if (sq == sr) then
          map%d(i,2) = nhelp1*nhelp2*(nhelp2-1)/2
        else
          map%d(i,2) = nhelp1*nhelp2*nhelp3
        end if
      else if ((typ == 7) .and. (sp == sq)) then
        map%d(i,2) = nhelp1*(nhelp1+1)*nhelp3/2
      else if ((typ == 8) .and. (sq == sr)) then
        map%d(i,2) = nhelp1*nhelp2*(nhelp2+1)/2
      else if (typ == 11) then
        if (sp == ss) then
          map%d(i,2) = nhelp1*(nhelp1+1)*(nhelp1+2)/6
        else if (sp == sq) then
          map%d(i,2) = nhelp1*(nhelp1+1)*nhelp3/2
        else if (sq == sr) then
          map%d(i,2) = nhelp1*nhelp2*(nhelp2+1)/2
        else
          map%d(i,2) = nhelp1*nhelp2*nhelp3
        end if
      else
        map%d(i,2) = nhelp1*nhelp2*nhelp3
      end if

      ! def sym p,q,r
      map%d(i,3) = sp
      map%d(i,4) = sq
      map%d(i,5) = sr
      map%d(i,6) = 0

      pos = pos+map%d(i,2)
      i = i+1

    end do
  end do

else if (nind == 4) then

  ! matrix A(p,q,r,s)

  do sp=1,nsym
    if ((typ == 1) .or. (typ == 4)) then
      nsymq = sp
    else
      nsymq = nsym
    end if

    do sq=1,nsymq
      spq = mul(sp,sq)
      if (typ == 2) then
        nsymr = sq
      else
        nsymr = nsym
      end if

      do sr=1,nsymr
        spqr = mul(spq,sr)

        ss = mul(stot,spqr)
        ! Meggie out
        if (((typ == 3) .or. (typ == 4)) .and. (sr < ss)) cycle

        nhelp1 = dimm(typp,sp)
        nhelp2 = dimm(typq,sq)
        nhelp3 = dimm(typr,sr)
        nhelp4 = dimm(typs,ss)

        ! def map%i
        map%i(sp,sq,sr) = i

        ! def position
        map%d(i,1) = pos

        ! def length
        if ((typ == 1) .and. (sp == sq)) then
          map%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
        else if ((typ == 2) .and. (sq == sr)) then
          map%d(i,2) = nhelp1*nhelp2*(nhelp3-1)*nhelp4/2
        else if ((typ == 3) .and. (sr == ss)) then
          map%d(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
        else if (typ == 4) then
          if ((sp == sq) .and. (sr == ss)) then
            map%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*(nhelp4-1)/4
          else if (sp == sq) then
            map%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
          else if (sr == ss) then
            map%d(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
          else
            map%d(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
          end if
        else
          map%d(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
        end if

        ! def sym p,q,r,s
        map%d(i,3) = sp
        map%d(i,4) = sq
        map%d(i,5) = sr
        map%d(i,6) = ss

        pos = pos+map%d(i,2)
        i = i+1

      end do
    end do
  end do

end if

post = pos

! definition of other coll

map%d(0,1) = typp
map%d(0,2) = typq
map%d(0,3) = typr
map%d(0,4) = typs
map%d(0,5) = i-1
map%d(0,6) = typ

return

end subroutine ccsort_t3grc0
