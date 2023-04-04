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

subroutine grc0stack(bsize,typ,typp,typq,typr,typs,stot,posst,map)
! This routine defines %d and %i for specific
! 3 index intermediate A(pq,Bp), needed when stacking
! (About Bp, see notes in multstack)
! This routine is a modification of grc0 routine
! P.N. 17.02.06
!
! N.B. (this routine cannot run with +OP2)

use ccsd_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: bsize, typ, typp, typq, typr, typs, stot, posst
type(Map_Type) :: map
integer(kind=iwp) :: i, nhelp1, nhelp2, nhelp3, poss, sp, sq

! To get rid of compiler warning
poss = 0
i = 0

! vanishing %i files

do nhelp1=1,nsym
  do nhelp2=1,nsym
    do nhelp3=1,nsym
      map%i(nhelp3,nhelp2,nhelp1) = 0
    end do
  end do
end do

! matrix A(p,q) or specifically A(i,j,Bp)

i = 1
poss = map%pos0

do sp=1,nsym

  sq = mmul(stot,sp)
  ! Meggie out
  if ((typ == 1) .and. (sp < sq)) cycle

  nhelp1 = dimm(typp,sp)
  nhelp2 = dimm(typq,sq)

  ! def %i
  map%i(sp,1,1) = i

  ! def position
  map%d(i,1) = poss

  ! def length
  if ((typ == 1) .and. (sp == sq)) then
    map%d(i,2) = bsize*nhelp1*(nhelp1-1)/2
  else
    map%d(i,2) = bsize*nhelp1*nhelp2
  end if

  ! def sym p,q
  map%d(i,3) = sp
  map%d(i,4) = sq
  map%d(i,5) = 0
  map%d(i,6) = 0

  poss = poss+map%d(i,2)
  i = i+1

end do

posst = poss

! definition of other coll

map%d(0,1) = typp
map%d(0,2) = typq
map%d(0,3) = typr
map%d(0,4) = typs
map%d(0,5) = i-1
map%d(0,6) = typ

return

end subroutine grc0stack
