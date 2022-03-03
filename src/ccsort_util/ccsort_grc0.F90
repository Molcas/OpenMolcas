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

subroutine ccsort_grc0(nind,typ,typp,typq,typr,typs,stot,pos0,post,mapd,mapi)
! this routine defines mapd and mapi for given intermediat

use ccsort_global, only: noa, nob, NSYM, nva, nvb
use Symmetry_Info, only: Mul
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nind, typ, typp, typq, typr, typs, stot, pos0
integer(kind=iwp), intent(out) :: post, mapd(0:512,6), mapi(8,8,8)
integer(kind=iwp) :: dimm(4,8), i, nhelp1, nhelp2, nhelp3, nhelp4, nsymq, nsymr, pos, sp, spq, spqr, sq, sr, ss

! !!!!!!!! def dimm to je tu len terazky !!!!

dimm(1,1:nsym) = noa(1:nsym)
dimm(2,1:nsym) = nob(1:nsym)
dimm(3,1:nsym) = nva(1:nsym)
dimm(4,1:nsym) = nvb(1:nsym)

! vanishing mapi files

mapi(1:nsym,1:nsym,1:nsym) = 0

!  To get rid of tiring compiler warning
i = 1
pos = pos0
if (nind == 1) then

  ! matrix A(p)
  sp = mul(stot,1)

  nhelp1 = dimm(typp,sp)

  ! def mapi
  mapi(1,1,1) = i

  ! def position
  mapd(i,1) = pos

  ! def length
  mapd(i,2) = nhelp1

  ! def sym p,q
  mapd(i,3) = sp
  mapd(i,4) = 0
  mapd(i,5) = 0
  mapd(i,6) = 0

  pos = pos+mapd(i,2)
  i = i+1

else if (nind == 2) then

  ! matrix A(p,q)

  do sp=1,nsym

    sq = mul(stot,sp)
    ! Meggie out
    if ((typ == 1) .and. (sp < sq)) cycle

    nhelp1 = dimm(typp,sp)
    nhelp2 = dimm(typq,sq)

    ! def mapi
    mapi(sp,1,1) = i

    ! def position
    mapd(i,1) = pos

    ! def length
    if ((typ == 1) .and. (sp == sq)) then
      mapd(i,2) = nhelp1*(nhelp1-1)/2
    else
      mapd(i,2) = nhelp1*nhelp2
    end if

    ! def sym p,q
    mapd(i,3) = sp
    mapd(i,4) = sq
    mapd(i,5) = 0
    mapd(i,6) = 0

    pos = pos+mapd(i,2)
    i = i+1

  end do

else if (nind == 3) then

  ! matrix A(p,q,r)

  do sp=1,nsym
    if (typ == 1) then
      nsymq = sp
    else
      nsymq = nsym
    end if

    do sq=1,nsymq
      spq = mul(sp,sq)

      sr = mul(stot,spq)
      ! Meggie out
      if ((typ == 2) .and. (sq < sr)) cycle

      nhelp1 = dimm(typp,sp)
      nhelp2 = dimm(typq,sq)
      nhelp3 = dimm(typr,sr)

      ! def mapi
      mapi(sp,sq,1) = i

      ! def position
      mapd(i,1) = pos

      ! def length
      if ((typ == 1) .and. (sp == sq)) then
        mapd(i,2) = nhelp1*(nhelp1-1)*nhelp3/2
      else if ((typ == 2) .and. (sq == sr)) then
        mapd(i,2) = nhelp1*nhelp2*(nhelp2-1)/2
      else
        mapd(i,2) = nhelp1*nhelp2*nhelp3
      end if

      ! def sym p,q,r
      mapd(i,3) = sp
      mapd(i,4) = sq
      mapd(i,5) = sr
      mapd(i,6) = 0

      pos = pos+mapd(i,2)
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

        ! def mapi
        mapi(sp,sq,sr) = i

        ! def position
        mapd(i,1) = pos

        ! def length
        if ((typ == 1) .and. (sp == sq)) then
          mapd(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
        else if ((typ == 2) .and. (sq == sr)) then
          mapd(i,2) = nhelp1*nhelp2*(nhelp3-1)*nhelp4/2
        else if ((typ == 3) .and. (sr == ss)) then
          mapd(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
        else if (typ == 4) then
          if ((sp == sq) .and. (sr == ss)) then
            mapd(i,2) = nhelp1*(nhelp2-1)*nhelp3*(nhelp4-1)/4
          else if (sp == sq) then
            mapd(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
          else if (sr == ss) then
            mapd(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
          else
            mapd(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
          end if
        else
          mapd(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
        end if

        ! def sym p,q,r,s
        mapd(i,3) = sp
        mapd(i,4) = sq
        mapd(i,5) = sr
        mapd(i,6) = ss

        pos = pos+mapd(i,2)
        i = i+1

      end do
    end do
  end do

end if

post = pos

! definition of other coll

mapd(0,1) = typp
mapd(0,2) = typq
mapd(0,3) = typr
mapd(0,4) = typs
mapd(0,5) = i-1
mapd(0,6) = typ

return

end subroutine ccsort_grc0
