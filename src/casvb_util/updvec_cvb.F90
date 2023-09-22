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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine updvec_cvb(upd,iorb,jorb,niprev,iprev,orbs,north,corth)
! Find update for IORB as projection of JORB on allowed space

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: upd(norb), orbs(norb,norb), corth(norb,niorth)
integer(kind=iwp) :: iorb, jorb, niprev, iprev(niprev), north(norb)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i1, io, ncon, noffort
real(kind=wp) :: dum(1)
integer(kind=iwp), external :: mstackr_cvb

i1 = mstackr_cvb(norb*norb)
noffort = 0
do io=1,iorb-1
  noffort = noffort+north(io)
end do
! Collect all constraints and find span:
call span0_cvb(norb,norb)
if (north(iorb) > 0) call span1_cvb(corth(1,1+noffort),north(iorb),dum,norb,0)
do i=1,niprev
  call span1_cvb(orbs(1,iprev(i)),1,dum,norb,0)
end do
call span1_cvb(orbs(1,iorb),1,dum,norb,0)
call span2_cvb(work(i1),ncon,dum,norb,0)

! Orthogonalise update to all remaining constraints
call fmove_cvb(orbs(1,jorb),upd,norb)
call schmidtd_cvb(work(i1),ncon,upd,1,dum,norb,0)
call mfreer_cvb(i1)

return

end subroutine updvec_cvb
