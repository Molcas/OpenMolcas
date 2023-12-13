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

subroutine cct3_t3grc0(nind,typ,typp,typq,typr,typs,stot,med,post)
! nind - number of indices (I)
! typ  - typ of mediate (I)
! typp - typ of index p (I)
! typq - typ of index q (I)
! typr - typ of index r (I)
! typs - typ of index s (I)
! stot - overall symmetry of the mediate (I)
! med  - mediate (I/O)
! post - final position of the mediate (O)
!
! this routine defines %d and %i for given intermediate
! it can done exactly the same maps like grc0 in CCSD
! plus additional types of mediates are introduced:
! type    meaning
! 5      p>q>r,s ; also p>q>r
! 6      p,q>r>s
! 7      p>=q,r,s ; also p>=q,r; p>=q
! 8      p,q>=r,s ; also p,q>=s
! 9      p,q,q>=s
! 10     p>=q,r>=s
! 11     p>=q>=r,s ; also p>=q>=r
! 12     p,q>=r>=s
!
! currently, these new types are implemented only for nind=3
!
! !N.B. (this routine cannot run with +OP2)
! N.B. this routine does not test stupidities

use CCT3_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nind, typ, typp, typq, typr, typs, stot
integer(kind=iwp), intent(out) :: post
type(Map_Type), intent(inout) :: med
integer(kind=iwp) :: i, nhelp1, nhelp2, nhelp3, nhelp4, nsymq, nsymr, pos, rsk1, rsk2, sp, spq, spqr, sq, sr, ss

! To fix some compiler warnings

ss = 0
pos = 0
rsk1 = 0
rsk2 = 0
i = 0

! vanishing med%i files

med%i(1:nsym,1:nsym,1:nsym) = 0

if (nind == 1) then

  ! matrix A(p)

  i = 1
  pos = med%pos0
  sp = mmul(stot,1)

  nhelp1 = dimm(typp,sp)

  ! def med%i
  med%i(1,1,1) = i

  ! def position
  med%d(i,1) = pos

  ! def length
  med%d(i,2) = nhelp1

  ! def sym p,q
  med%d(i,3) = sp
  med%d(i,4) = 0
  med%d(i,5) = 0
  med%d(i,6) = 0

  pos = pos+med%d(i,2)
  i = i+1

else if (nind == 2) then

  ! matrix A(p,q)

  i = 1
  pos = med%pos0

  do sp=1,nsym

    sq = mmul(stot,sp)
    ! Meggie out
    if ((typ == 1) .and. (sp < sq)) cycle

    nhelp1 = dimm(typp,sp)
    nhelp2 = dimm(typq,sq)

    ! def med%i
    med%i(sp,1,1) = i

    ! def position
    med%d(i,1) = pos

    ! def length
    if ((typ == 1) .and. (sp == sq)) then
      med%d(i,2) = nhelp1*(nhelp1-1)/2
    else
      med%d(i,2) = nhelp1*nhelp2
    end if

    ! def sym p,q
    med%d(i,3) = sp
    med%d(i,4) = sq
    med%d(i,5) = 0
    med%d(i,6) = 0

    pos = pos+med%d(i,2)
    i = i+1

  end do

else if (nind == 3) then

  ! matrix A(p,q,r)

  ! def reucion sumations keys : rsk1 for pq, rsk2 for qr

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

  i = 1
  pos = med%pos0

  do sp=1,nsym
    if (rsk1 == 1) then
      nsymq = sp
    else
      nsymq = nsym
    end if

    do sq=1,nsymq
      spq = mmul(sp,sq)

      sr = mmul(stot,spq)
      ! Meggie out
      if ((rsk2 == 1) .and. (sq < sr)) cycle

      nhelp1 = dimm(typp,sp)
      nhelp2 = dimm(typq,sq)
      nhelp3 = dimm(typr,sr)

      ! def med%i
      med%i(sp,sq,1) = i

      ! def position
      med%d(i,1) = pos

      ! def length
      if ((typ == 1) .and. (sp == sq)) then
        med%d(i,2) = nhelp1*(nhelp1-1)*nhelp3/2
      else if ((typ == 2) .and. (sq == sr)) then
        med%d(i,2) = nhelp1*nhelp2*(nhelp2-1)/2
      else if (typ == 5) then
        if (sp == sr) then
          med%d(i,2) = nhelp1*(nhelp1-1)*(nhelp1-2)/6
        else if (sp == sq) then
          med%d(i,2) = nhelp1*(nhelp1-1)*nhelp3/2
        else if (sq == sr) then
          med%d(i,2) = nhelp1*nhelp2*(nhelp2-1)/2
        else
          med%d(i,2) = nhelp1*nhelp2*nhelp3
        end if
      else if ((typ == 7) .and. (sp == sq)) then
        med%d(i,2) = nhelp1*(nhelp1+1)*nhelp3/2
      else if ((typ == 8) .and. (sq == sr)) then
        med%d(i,2) = nhelp1*nhelp2*(nhelp2+1)/2
      else if (typ == 11) then
        if (sp == ss) then
          med%d(i,2) = nhelp1*(nhelp1+1)*(nhelp1+2)/6
        else if (sp == sq) then
          med%d(i,2) = nhelp1*(nhelp1+1)*nhelp3/2
        else if (sq == sr) then
          med%d(i,2) = nhelp1*nhelp2*(nhelp2+1)/2
        else
          med%d(i,2) = nhelp1*nhelp2*nhelp3
        end if
      else
        med%d(i,2) = nhelp1*nhelp2*nhelp3
      end if

      ! def sym p,q,r
      med%d(i,3) = sp
      med%d(i,4) = sq
      med%d(i,5) = sr
      med%d(i,6) = 0

      pos = pos+med%d(i,2)
      i = i+1

    end do
  end do

else if (nind == 4) then

  ! matrix A(p,q,r,s)

  i = 1
  pos = med%pos0

  do sp=1,nsym
    if ((typ == 1) .or. (typ == 4)) then
      nsymq = sp
    else
      nsymq = nsym
    end if

    do sq=1,nsymq
      spq = mmul(sp,sq)
      if (typ == 2) then
        nsymr = sq
      else
        nsymr = nsym
      end if

      do sr=1,nsymr
        spqr = mmul(spq,sr)

        ss = mmul(stot,spqr)
        ! Meggie out
        if (((typ == 3) .or. (typ == 4)) .and. (sr < ss)) cycle

        nhelp1 = dimm(typp,sp)
        nhelp2 = dimm(typq,sq)
        nhelp3 = dimm(typr,sr)
        nhelp4 = dimm(typs,ss)

        ! def med%i
        med%i(sp,sq,sr) = i

        ! def position
        med%d(i,1) = pos

        ! def length
        if ((typ == 1) .and. (sp == sq)) then
          med%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
        else if ((typ == 2) .and. (sq == sr)) then
          med%d(i,2) = nhelp1*nhelp2*(nhelp3-1)*nhelp4/2
        else if ((typ == 3) .and. (sr == ss)) then
          med%d(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
        else if (typ == 4) then
          if ((sp == sq) .and. (sr == ss)) then
            med%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*(nhelp4-1)/4
          else if (sp == sq) then
            med%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
          else if (sr == ss) then
            med%d(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
          else
            med%d(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
          end if
        else
          med%d(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
        end if

        ! def sym p,q,r,s
        med%d(i,3) = sp
        med%d(i,4) = sq
        med%d(i,5) = sr
        med%d(i,6) = ss

        pos = pos+med%d(i,2)
        i = i+1

      end do
    end do
  end do

end if

post = pos

! definition of other coll

med%d(0,1) = typp
med%d(0,2) = typq
med%d(0,3) = typr
med%d(0,4) = typs
med%d(0,5) = i-1
med%d(0,6) = typ

return

end subroutine cct3_t3grc0
