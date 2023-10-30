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

subroutine serber_cvb(bikcof,nel,nalf,nbet,ndet,ifns)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nel, nalf, nbet, ndet, ifns
real(kind=wp), intent(inout) :: bikcof(ndet,ifns)
integer(kind=iwp) :: ia, iachek, iantival, ib, indx, iorb, k, l, nc, rc
integer(kind=iwp), allocatable :: ialfs(:), ianti(:), ibets(:), lnocca(:), locca(:), maxspn(:), minspn(:), nkspn(:), xspin(:)
real(kind=wp) :: dum(1)

! Serber spin functions

! For each Rumer spin function we determine (minus) the "antisymmetry"
! number, as explained in [Spin Eigenfunctions, R. Pauncz, Sec. 5.5].

call mma_allocate(minspn,[0,nel],label='minspn')
call mma_allocate(maxspn,[0,nel],label='maxspn')
call mma_allocate(nkspn,[0,nel],label='nkspn')
call mma_allocate(locca,nel,label='locca')
call mma_allocate(lnocca,nel,label='lnocca')
call mma_allocate(xspin,(nel+1)*(nalf+1),label='xspin')
call mma_allocate(ialfs,nalf,label='ialfs')
call mma_allocate(ibets,nbet,label='ibets')
call mma_allocate(ianti,ifns,label='ianti')

! Spin function weight arrays:
do iorb=0,nel
  minspn(iorb) = max(iorb-nalf,0)
  maxspn(iorb) = min(iorb/2,nbet)
end do
call weight_cvb(xspin,minspn,maxspn,nbet,nel)
if (ifns /= xspin((nel+1)*(nbet+1))) then
  write(u6,*) ' Discrepancy in IFNS:',ifns,xspin((nel+1)*(nbet+1))
  call abend_cvb()
end if
nkspn(:) = maxspn(:)
call occupy_cvb(nkspn,nel,locca,lnocca)
! Loop:
indx = 1
! Determine pairings
do
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

  ianti(indx) = 0
  do ib=1,nbet
    if ((mod(ialfs(ib),2) == 1) .and. (ialfs(ib) == ibets(ib)-1)) ianti(indx) = ianti(indx)-1
  end do
  call loind_cvb(nel,nbet,nkspn,minspn,maxspn,locca,lnocca,indx,xspin,rc)
  if (rc == 0) exit
end do

call mma_deallocate(minspn)
call mma_deallocate(maxspn)
call mma_deallocate(nkspn)
call mma_deallocate(locca)
call mma_deallocate(lnocca)
call mma_deallocate(xspin)
call mma_deallocate(ialfs)
call mma_deallocate(ibets)

! Sort according to decreasing values of IANTI:
nc = 0
do iantival=nbet,0,-1
  do k=1,ifns
    if (ianti(k) == -iantival) then
      nc = nc+1
      ianti(k) = nc
    end if
  end do
end do

do k=1,ifns
  if (ianti(k) /= k) then
    do l=1,ifns
      if (ianti(l) == k) exit
    end do
    if (l > ifns) then
      write(u6,*) ' Error - swap function not found!',k,ianti(k)
      call abend_cvb()
    end if
    call dswap_(ndet,bikcof(:,k),1,bikcof(:,l),1)
    ianti(l) = ianti(k)
    ianti(k) = k
  end if
end do

call mma_deallocate(ianti)

! Orthonormalize:
call schmidtn_cvb(bikcof,ifns,dum,ndet,0)

return

end subroutine serber_cvb
