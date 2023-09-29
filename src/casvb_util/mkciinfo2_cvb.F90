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

subroutine mkciinfo2_cvb(i1alf,i1bet,iafrm,ibfrm,iato,ibto,phato,phbto)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: i1alf(n1a,norb), i1bet(n1b,norb), iafrm(norb,nda), ibfrm(norb,ndb), iato(norb,0:nam1), ibto(norb,0:nbm1)
real(kind=wp) :: phato(norb,nam1), phbto(norb,nbm1)
integer(kind=iwp) :: i, ia, iax, iaxtmp, ib, ibx, ibxtmp, iel, indx, iorb, rc
integer(kind=iwp), allocatable :: iaccm(:), inewocc(:), locc(:), lunocc(:), maxgrph(:), mingrph(:), nk(:), xalf(:,:), xalf2(:,:), &
                                  xbet(:,:), xbet2(:,:)
integer(kind=iwp), external :: indget_cvb, ip_of_iWork

call mma_allocate(xalf,[0,norb],[0,nalf],label='xalf')
call mma_allocate(xbet,[0,norb],[0,nbet],label='xbet')
call mma_allocate(xalf2,[0,norb],[0,nalf-1],label='xalf2')
call mma_allocate(xbet2,[0,norb],[0,nbet-1],label='xbet2')
call mma_allocate(mingrph,[0,norb],label='mingrph')
call mma_allocate(maxgrph,[0,norb],label='maxgrph')
call mma_allocate(nk,[0,norb],label='nk')
call mma_allocate(locc,norb+1,label='locc')
call mma_allocate(lunocc,norb+1,label='lunocc')
call mma_allocate(inewocc,norb,label='inewocc')
call mma_allocate(iaccm,norb,label='iaccm')

call izero(iafrm,nda*norb)
call izero(ibfrm,ndb*norb)
call izero(iato,(nam1+1)*norb)
call izero(ibto,(nbm1+1)*norb)
call fzero(phato,nam1*norb)
call fzero(phbto,nbm1*norb)
! Alpha loop:
call izero(iaccm,norb)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf,0)
  maxgrph(iorb) = min(iorb,nalf)
end do
call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
call imove_cvb(maxgrph,nk,norb+1)
call occupy_cvb(nk,norb,locc,lunocc)
indx = 1
do
  call izero(inewocc,norb)
  do i=1,nalf
    inewocc(locc(i)) = 1
  end do
  do iel=1,norb
    if (inewocc(iel) == 1) then
      iaccm(iel) = iaccm(iel)+1
      i1alf(iaccm(iel),iel) = indx
    end if
  end do
  ! Indexing arrays:
  do i=1,nalf
    inewocc(locc(i)) = 0
    iafrm(locc(i),indx) = indget_cvb(inewocc,nalf,norb,xalf)
    inewocc(locc(i)) = 1
  end do
  call loind_cvb(norb,nalf,nk,mingrph,maxgrph,locc,lunocc,indx,xalf,rc)
  if (rc == 0) exit
end do
! Beta loop:
call izero(iaccm,norb)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nbet,0)
  maxgrph(iorb) = min(iorb,nbet)
end do
call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)
call imove_cvb(maxgrph,nk,norb+1)
call occupy_cvb(nk,norb,locc,lunocc)
indx = 1
do
  call izero(inewocc,norb)
  do i=1,nbet
    inewocc(locc(i)) = 1
  end do
  do iel=1,norb
    if (inewocc(iel) == 1) then
      iaccm(iel) = iaccm(iel)+1
      i1bet(iaccm(iel),iel) = indx
    end if
  end do
  ! Indexing arrays:
  do i=1,nbet
    inewocc(locc(i)) = 0
    ibfrm(locc(i),indx) = indget_cvb(inewocc,nbet,norb,xbet)
    inewocc(locc(i)) = 1
  end do
  call loind_cvb(norb,nbet,nk,mingrph,maxgrph,locc,lunocc,indx,xbet,rc)
  if (rc == 0) exit
end do

! Altered definitions of I1ALF & I1BET:
do iorb=1,norb
  do ia=1,n1a
    iax = i1alf(ia,iorb)
    iaxtmp = iafrm(iorb,iax)
    i1alf(ia,iorb) = iaxtmp
  end do
end do
if (absym(4)) then
  ! I1ALF & I1BET may share memory:
  if (ip_of_iWork(i1alf(1,1)) /= ip_of_iWork(i1bet(1,1))) call imove_cvb(i1alf,i1bet,norb*n1a)
else
  do iorb=1,norb
    do ib=1,n1b
      ibx = i1bet(ib,iorb)
      ibxtmp = ibfrm(iorb,ibx)
      i1bet(ib,iorb) = ibxtmp
    end do
  end do
end if
! More indexing arrays:
! Second alpha loop (NALF-1):
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf-1,0)
  maxgrph(iorb) = min(iorb,nalf-1)
end do
call weight_cvb(xalf2,mingrph,maxgrph,nalf-1,norb)
call imove_cvb(maxgrph,nk,norb+1)
call occupy_cvb(nk,norb,locc,lunocc)
indx = 1
do
  call izero(inewocc,norb)
  do i=1,nalf-1
    inewocc(locc(i)) = 1
  end do
  do i=1,norb-nalf+1
    inewocc(lunocc(i)) = 1
    iato(lunocc(i),indx) = indget_cvb(inewocc,nalf,norb,xalf)
    inewocc(lunocc(i)) = 0
    if (mod(nalf+i+1+lunocc(i),2) == 0) then
      phato(lunocc(i),indx) = one
    else
      phato(lunocc(i),indx) = -one
    end if
  end do
  call loind_cvb(norb,nalf-1,nk,mingrph,maxgrph,locc,lunocc,indx,xalf2,rc)
  if (rc == 0) exit
end do

if (nbet > 0) then
  ! Second beta loop (NBET-1):
  do iorb=0,norb
    mingrph(iorb) = max(iorb-norb+nbet-1,0)
    maxgrph(iorb) = min(iorb,nbet-1)
  end do
  call weight_cvb(xbet2,mingrph,maxgrph,nbet-1,norb)
  call imove_cvb(maxgrph,nk,norb+1)
  call occupy_cvb(nk,norb,locc,lunocc)
  indx = 1
  do
    call izero(inewocc,norb)
    do i=1,nbet-1
      inewocc(locc(i)) = 1
    end do
    do i=1,norb-nbet+1
      inewocc(lunocc(i)) = 1
      ibto(lunocc(i),indx) = indget_cvb(inewocc,nbet,norb,xbet)
      inewocc(lunocc(i)) = 0
      if (mod(nbet+i+1+lunocc(i),2) == 0) then
        phbto(lunocc(i),indx) = one
      else
        phbto(lunocc(i),indx) = -one
      end if
    end do
    call loind_cvb(norb,nbet-1,nk,mingrph,maxgrph,locc,lunocc,indx,xbet2,rc)
    if (rc == 0) exit
  end do
end if

call mma_deallocate(xalf)
call mma_deallocate(xbet)
call mma_deallocate(xalf2)
call mma_deallocate(xbet2)
call mma_deallocate(mingrph)
call mma_deallocate(maxgrph)
call mma_deallocate(nk)
call mma_deallocate(locc)
call mma_deallocate(lunocc)
call mma_deallocate(inewocc)
call mma_deallocate(iaccm)

return

end subroutine mkciinfo2_cvb
