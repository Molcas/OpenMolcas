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

subroutine mkorbfree_cvb()

use casvb_global, only: corth, ifxorb, iorts, irels, irots, nfr, nfrorb, nfrvb, nijrel, norb, north, nprorb, orbfr_is_unit, &
                        orbopt, orbs, relorb, trprm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, i2, icon, ierr, ioff, ioff2, ioffs, iorb, irel, ishift, ishift2, j, j2, jorb, nc, ncon, nl1, nl2, nrem
real(kind=wp) :: dum(1), sum1, sum2
integer(kind=iwp), allocatable :: idel(:)
real(kind=wp), allocatable :: orbinv(:,:), owrk(:,:), owrk2(:,:)
real(kind=wp), parameter :: thresh = 1.0e-7_wp
real(kind=wp), external :: ddot_

if (orbfr_is_unit) then
  nfrorb = nprorb
else

  call mma_allocate(owrk,norb,norb,label='owrk')
  call mma_allocate(owrk2,norb,norb,label='owrk2')
  call mma_allocate(orbinv,norb,norb,label='orbinv')
  call mma_allocate(idel,nprorb,label='idel')

  trprm(:,:) = Zero
  idel(:) = 0
  orbinv(:,:) = orbs(:,:)
  call mxinv_cvb(orbinv,norb)

  nc = 0
  ioffs = 0
  do iorb=1,norb
    if ((north(iorb) > 0) .and. (ifxorb(iorb) /= 1)) then
      ! Transform simple constraints to basis of VB orbitals:
      call mxattb_cvb(orbs,corth(:,ioffs+1:),norb,norb,north(iorb),owrk)
      call span_cvb(owrk,north(iorb),ncon,dum,norb,0)
      ishift = (iorb-1)*(norb-1)
      do icon=1,ncon
        nc = nc+1
        i2 = 0
        do i=1,norb
          if (i /= iorb) then
            i2 = i2+1
            trprm(i2+ishift,nc) = owrk(i,icon)
          end if
        end do
      end do
    else if (ifxorb(iorb) == 1) then
      ishift = (iorb-1)*(norb-1)
      do icon=1,norb-1
        nc = nc+1
        trprm(icon+ishift,nc) = One
      end do
    end if
    ioffs = ioffs+north(iorb)
  end do

  call mxattb_cvb(orbs,orbs,norb,norb,norb,owrk)
  call ortelim_cvb(trprm,iorts,irots,owrk,nc,nprorb,norb*(norb-1),nrem)
  idel(1:nrem) = 1
  idel(nrem+1:) = 0

  do irel=1,nijrel
    iorb = irels(1,irel)
    jorb = irels(2,irel)
    call mxatb_cvb(relorb(:,:,irel),orbs,norb,norb,norb,owrk)
    call mxatb_cvb(orbinv,owrk,norb,norb,norb,owrk2)
    if (abs(abs(owrk2(iorb,jorb))-One) > thresh) then
      write(u6,*) ' Transformation matrix cannot be correct !'
      call mxprint_cvb(owrk2,norb,norb,0)
      call abend_cvb()
    end if
    ishift = (jorb-1)*(norb-1)
    ishift2 = (iorb-1)*(norb-1)
    i2 = 0
    do i=1,norb
      if (i == iorb) cycle
      i2 = i2+1
      j2 = 0
      do j=1,norb
        if (j == jorb) cycle
        j2 = j2+1
        trprm(i2+ishift2,:) = trprm(i2+ishift2,:)+owrk2(i,j)*trprm(j2+ishift,:)
      end do
    end do
    ioff = 1+(iorb-1)*(norb-1)
    ioff2 = 1+iorb*(norb-1)
    nl1 = (iorb-1)*(norb-1)
    nl2 = (norb-iorb)*(norb-1)
    do i=nrem+1,nprorb
      sum1 = ddot_(norb-1,trprm(ioff:,i),1,trprm(ioff:,i),1)
      sum2 = ddot_(nl1,trprm(:,i),1,trprm(:,i),1)
      if (nl2 > 0) sum2 = sum2+ddot_(nl2,trprm(ioff2:,i),1,trprm(ioff2:,i),1)
      if ((sum1 > thresh) .and. (sum2 < thresh)) idel(i) = 1
    end do
  end do
  nfrorb = 0
  do i=1,norb*(norb-1)
    if (idel(i) /= 1) then
      nfrorb = nfrorb+1
      trprm(:,nfrorb) = trprm(:,i)
    end if
  end do
  trprm(:,nfrorb+1:) = Zero
  ierr = 0
  call nize_cvb(trprm,nfrorb,dum,nprorb,0,ierr)

  call mma_deallocate(owrk)
  call mma_deallocate(owrk2)
  call mma_deallocate(orbinv)
  call mma_deallocate(idel)

end if

nfr = nfrorb+nfrvb
orbopt = (nfrorb /= 0)

return

end subroutine mkorbfree_cvb
