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

subroutine ortelim_cvb(trprm,iorts,irots,sorbs,nc,npr1,norbprm,nrem)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: iorts(2,nort), irots(2,ndrot), nc, npr1, norbprm, nrem
real(kind=wp) :: trprm(npr1,npr1), sorbs(norb,norb)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i1, i1ff, iorb, iort, irot, jorb, ki, kj, korb
real(kind=wp) :: dum(1)
integer(kind=iwp), external :: mstackrz_cvb

i1 = mstackrz_cvb(norbprm*max(nc+nort+ndrot,norbprm))
i1ff = i1-1
do i=1,nc
  call fmove_cvb(trprm(1,i),work(1+(i-1)*norbprm+i1ff),norbprm)
end do
do iort=1,nort
  iorb = iorts(1,iort)
  jorb = iorts(2,iort)
  do korb=1,norb
    ki = korb+(iorb-1)*(norb-1)
    if (korb > iorb) ki = ki-1
    kj = korb+(jorb-1)*(norb-1)
    if (korb > jorb) kj = kj-1
    if (korb /= iorb) work(ki+(iort+nc-1)*norbprm+i1ff) = sorbs(korb,jorb)
    if (korb /= jorb) work(kj+(iort+nc-1)*norbprm+i1ff) = sorbs(korb,iorb)
  end do
end do
do irot=1,ndrot
  iorb = irots(1,irot)
  jorb = irots(2,irot)
  do korb=1,norb
    ki = korb+(iorb-1)*(norb-1)
    if (korb > iorb) ki = ki-1
    kj = korb+(jorb-1)*(norb-1)
    if (korb > jorb) kj = kj-1
    if (korb /= iorb) work(ki+(irot+nc+nort-1)*norbprm+i1ff) = sorbs(korb,jorb)
    if (korb /= jorb) work(kj+(irot+nc+nort-1)*norbprm+i1ff) = -sorbs(korb,iorb)
  end do
end do
call span_cvb(work(i1),nc+nort+ndrot,nrem,dum,norbprm,0)
call compl_cvb(work(i1),nrem,norbprm)

call fzero(trprm,npr1*npr1)
do i=1,norbprm
  call fmove_cvb(work(1+(i-1)*norbprm+i1ff),trprm(1,i),norbprm)
end do

call mfreer_cvb(i1)

return

end subroutine ortelim_cvb
