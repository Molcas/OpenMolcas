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

subroutine setipermzeta_cvb(ipermzeta,orbs,symelm,izeta)

use casvb_global, only: norb, nsyme, nzeta, tags
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ipermzeta(norb,nzeta)
real(kind=wp), intent(in) :: orbs(norb,norb), symelm(norb*norb,nsyme)
integer(kind=iwp), intent(in) :: izeta(nsyme)
integer(kind=iwp) :: iorb, isyme, izeta1, jorb
real(kind=wp), allocatable :: orbinv(:,:), owrk(:,:), owrk2(:,:)
real(kind=wp), parameter :: thresh = 1.0e-8_wp

call mma_allocate(orbinv,norb,norb,label='orbinv')
call mma_allocate(owrk,norb,norb,label='owrk')
call mma_allocate(owrk2,norb,norb,label='owrk2')

if (nzeta > 0) then
  orbinv(:,:) = orbs(:,:)
  call mxinv_cvb(orbinv,norb)
end if

izeta1 = 0
do isyme=1,nsyme
  if (izeta(isyme) /= 0) then
    izeta1 = izeta1+1
    ! Determine orbital permutation for sym. operation ISYME:
    call mxatb_cvb(symelm(:,isyme),orbs,norb,norb,norb,owrk2)
    call mxatb_cvb(orbinv,owrk2,norb,norb,norb,owrk)
    do iorb=1,norb
      do jorb=1,norb
        if (abs(abs(owrk(jorb,iorb))-One) < thresh) then
          ipermzeta(iorb,izeta1) = nint(owrk(jorb,iorb))*jorb
        else if (abs(owrk(jorb,iorb)) > thresh) then
          write(u6,*) ' Fatal error! Symmetry operation ',tags(isyme),' does not permute the VB orbitals!'
          call mxprint_cvb(owrk,norb,norb,0)
          call abend_cvb()
        end if
      end do
    end do
  end if
end do

call mma_deallocate(orbinv)
call mma_deallocate(owrk)
call mma_deallocate(owrk2)

return

end subroutine setipermzeta_cvb
