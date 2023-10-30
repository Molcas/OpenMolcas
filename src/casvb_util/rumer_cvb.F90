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

subroutine rumer_cvb(bikcof,nel,nalf,nbet,ndet,ifns,kbasis,iprint,nswpdim)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nel, nalf, nbet, ndet, ifns, kbasis, iprint, nswpdim
real(kind=wp), intent(out) :: bikcof(ndet,ifns)
integer(kind=iwp) :: ia, iachek, ib, ii, indx, iorb, iswp, nbet2, rc
real(kind=wp) :: abphase, bikvalue, scl
integer(kind=iwp), allocatable :: ialfs(:), ibets(:), ioccswp(:,:), iw(:), locca(:), lnocca(:), maxspn(:), maxswp(:), minspn(:), &
                                  minswp(:), nkspn(:), nkswp(:), xdet(:,:), xspin(:)
integer(kind=iwp), external :: indget_cvb
real(kind=wp), external :: party_cvb

bikcof(:,:) = Zero

call mma_allocate(xdet,[0,nel],[0,nalf],label='xdet')

! Determinant weight array (index from alpha spin string):
call weightfl_cvb(xdet,nalf,nel)
if (ndet /= xdet(nel,nalf)) then
  write(u6,*) ' Discrepancy in NDET:',ndet,xdet(nel,nalf)
  call abend_cvb()
end if

! Rumer spin functions

if ((iprint >= 2) .and. ((kbasis == 3) .or. (kbasis == 4))) write(u6,6100) 2**nbet
abphase = real(1-2*mod(nTri_Elem(nbet-1),2),kind=wp)

call mma_allocate(minswp,[0,2*nbet],label='minswp')
call mma_allocate(maxswp,[0,2*nbet],label='maxswp')
call mma_allocate(nkswp,[0,2*nbet],label='nkswp')
call mma_allocate(ioccswp,nbet,nswpdim,label='ioccswp')

! Prepare NKs for a<->b interchanges in (ab-ba) terms:
nbet2 = 2*nbet
do iorb=0,nbet2
  minswp(iorb) = iorb/2
  maxswp(iorb) = (iorb+1)/2
end do
nkswp(:) = maxswp(:)
iswp = 0
do
  iswp = iswp+1
  do ia=1,nbet
    ioccswp(ia,iswp) = nkswp(2*ia)-nkswp(2*ia-1)
  end do
  call loop_cvb(nbet2,nkswp,minswp,maxswp,rc)
  if (rc == 0) exit
end do

call mma_deallocate(minswp)
call mma_deallocate(maxswp)
call mma_deallocate(nkswp)
call mma_allocate(minspn,[0,nel],label='minspn')
call mma_allocate(maxspn,[0,nel],label='maxspn')
call mma_allocate(nkspn,[0,nel],label='nkspn')
call mma_allocate(locca,nel,label='locca')
call mma_allocate(lnocca,nel,label='lnocca')
call mma_allocate(xspin,(nel+1)*(nalf+1),label='xspin')

! Spin function weight arrays:
do iorb=0,nel
  minspn(iorb) = max(iorb-nalf,0)
  maxspn(iorb) = min(iorb/2,nbet)
end do
call weight_cvb(xspin,minspn,maxspn,nbet,nel)
if ((ifns /= xspin((nel+1)*(nbet+1))) .and. (kbasis /= 6)) then
  write(u6,*) ' Discrepancy in IFNS:',ifns,xspin((nel+1)*(nbet+1))
  call abend_cvb()
end if
nkspn(:) = maxspn(:)
call occupy_cvb(nkspn,nel,locca,lnocca)

call mma_allocate(iw,nel,label='iw')
call mma_allocate(ialfs,nalf,label='ialfs')
call mma_allocate(ibets,nbet,label='ibets')

! Loop:
indx = 1
! Determine pairings
do while ((kbasis /= 6) .or. (indx <= ifns))
  do ib=1,nbet
    ibets(ib) = locca(ib)
    do ia=nalf,1,-1
      ialfs(ib) = lnocca(ia)
      if (ialfs(ib) < ibets(ib)) then
        do iachek=1,ib-1
          if (ialfs(iachek) == ialfs(ib)) exit
        end do
        if (iachek >= ib) exit
      end if
    end do
  end do
  do ib=nbet+1,nalf
    do ia=1,nalf
      ialfs(ib) = lnocca(ia)
      do iachek=1,ib-1
        if (ialfs(iachek) == ialfs(ib)) exit
      end do
      if (iachek >= ib) exit
    end do
  end do

  if ((iprint >= 2) .and. ((kbasis == 3) .or. (kbasis == 4))) write(u6,6200) indx,(ialfs(ii),ibets(ii),ii=1,nbet)

  if (kbasis == 4) then
    bikvalue = One
  else
    bikvalue = abphase*party_cvb(ialfs,nalf)*party_cvb(ibets,nbet)
  end if

  do iswp=1,nswpdim
    iw(:) = 0
    do ia=1,nalf
      iw(ialfs(ia)) = 1
    end do
    do ia=1,nbet
      if (ioccswp(ia,iswp) /= 0) then
        ! IA => alpha electron in position number 2 (= swap).
        iw(ialfs(ia)) = 0
        iw(ibets(ia)) = 1
      end if
    end do
    bikcof(indget_cvb(iw,nalf,nel,xdet),indx) = bikvalue
  end do
  call loind_cvb(nel,nbet,nkspn,minspn,maxspn,locca,lnocca,indx,xspin,rc)
  if (rc == 0) exit
end do

call mma_deallocate(minspn)
call mma_deallocate(maxspn)
call mma_deallocate(nkspn)
call mma_deallocate(ioccswp)
call mma_deallocate(locca)
call mma_deallocate(lnocca)
call mma_deallocate(xdet)
call mma_deallocate(xspin)
call mma_deallocate(iw)
call mma_deallocate(ialfs)
call mma_deallocate(ibets)

! Normalise
scl = One/sqrt(real(2**nbet,kind=wp))
bikcof(:,:) = scl*bikcof(:,:)

return

6100 format(/,' Number of determinants per structure:',i4)
6200 format(2x,i3,' ==> ',8(i2,'-',i2,3x))

end subroutine rumer_cvb
