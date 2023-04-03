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

subroutine grc0stack(bsize,typ,typp,typq,typr,typs,stot,poss0,posst,mapd,mapi)
! This routine defines mapd and mapi for specific
! 3 index intermediate A(pq,Bp), needed when stacking
! (About Bp, see notes in multstack)
! This routine is a modification of grc0 routine
! P.N. 17.02.06
!
! N.B. (this routine cannot run with +OP2)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: bsize, typ, typp, typq, typr, typs, stot, poss0, posst, mapd(0:512,6), mapi(8,8,8)
#include "ccsd1.fh"
integer(kind=iwp) :: i, nhelp1, nhelp2, nhelp3, poss, sp, sq

! To get rid of compiler warning
poss = 0
i = 0

! vanishing mapi files

do nhelp1=1,nsym
  do nhelp2=1,nsym
    do nhelp3=1,nsym
      mapi(nhelp3,nhelp2,nhelp1) = 0
    end do
  end do
end do

! matrix A(p,q) or specifically A(i,j,Bp)

i = 1
poss = poss0

do sp=1,nsym

  sq = mmul(stot,sp)
  ! Meggie out
  if ((typ == 1) .and. (sp < sq)) cycle

  nhelp1 = dimm(typp,sp)
  nhelp2 = dimm(typq,sq)

  ! def mapi
  mapi(sp,1,1) = i

  ! def position
  mapd(i,1) = poss

  ! def length
  if ((typ == 1) .and. (sp == sq)) then
    mapd(i,2) = bsize*nhelp1*(nhelp1-1)/2
  else
    mapd(i,2) = bsize*nhelp1*nhelp2
  end if

  ! def sym p,q
  mapd(i,3) = sp
  mapd(i,4) = sq
  mapd(i,5) = 0
  mapd(i,6) = 0

  poss = poss+mapd(i,2)
  i = i+1

end do

posst = poss

! definition of other coll

mapd(0,1) = typp
mapd(0,2) = typq
mapd(0,3) = typr
mapd(0,4) = typs
mapd(0,5) = i-1
mapd(0,6) = typ

return

end subroutine grc0stack
