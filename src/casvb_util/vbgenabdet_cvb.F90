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

subroutine vbgenabdet_cvb(idetavb,idetbvb,iconfs,nconf,nconfion,ndetvb,nel,noe,nalf,nbet,norb)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nconf, noe, iconfs(noe,nconf), nel, nconfion(0:nel), ndetvb, nalf, nbet, norb
integer(kind=iwp), intent(out) :: idetavb(ndetvb), idetbvb(ndetvb)
integer(kind=iwp) :: i, iaind, iaocc, ibind, ibocc, iconf, incr, incrdet, indx, ioff_nconf, ion, iorb, nalfsing, nbetsing, &
                     nelsing, nstring
integer(kind=iwp), allocatable :: astr(:,:), bstr(:,:), iaccm(:), inewocc(:), maxgrph(:), mingrph(:), xalf(:,:), xbet(:,:)
logical(kind=iwp), parameter :: debug = .false.
integer(kind=iwp), external :: indget_cvb

call mma_allocate(xalf,[0,norb],[0,nalf],label='xalf')
call mma_allocate(xbet,[0,norb],[0,nbet],label='xbet')
call mma_allocate(mingrph,[0,norb],label='mingrph')
call mma_allocate(maxgrph,[0,norb],label='maxgrph')

! Set xalf and xbet for indget:
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf,0)
  maxgrph(iorb) = min(iorb,nalf)
end do
call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nbet,0)
  maxgrph(iorb) = min(iorb,nbet)
end do
call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)

call mma_deallocate(mingrph)
call mma_deallocate(maxgrph)
call mma_allocate(inewocc,norb,label='inewocc')
call mma_allocate(iaccm,norb,label='iaccm')

! Transformation matrix VB structures  <->  determinants
incrdet = 0
ioff_nconf = 0
do ion=0,nel/2
  nelsing = nel-2*ion
  nalfsing = nalf-ion
  nbetsing = nbet-ion
  ! Skip if configurations are incompatible with Ms:
  if ((nalfsing >= 0) .and. (nbetsing >= 0) .and. (nelsing >= 0)) then

    ! Generate alpha/beta strings for singly occupied electrons:
    call icomb_cvb(nelsing,nalfsing,nstring)
    call mma_allocate(astr,nalfsing,nstring,label='astr')
    call mma_allocate(bstr,nbetsing,nstring,label='bstr')
    call stringen_cvb(nelsing,nalfsing,astr,bstr)
    if (debug) then
      write(u6,*) ' ionicity=',ion,' nconf=',nconfion(ion)
      write(u6,*) ' check alpha strings :'
      do i=1,nstring
        write(u6,*) i,' => ',astr(:,i)
      end do
      write(u6,*) ' check beta strings :'
      do i=1,nstring
        write(u6,*) i,' => ',bstr(:,i)
      end do
    end if

    do iconf=ioff_nconf+1,ioff_nconf+nconfion(ion)
      inewocc(:) = iconfs(1:norb,iconf)
      incr = 0
      do iorb=1,norb
        if (inewocc(iorb) == 1) then
          incr = incr+1
          iaccm(incr) = iorb
        end if
        inewocc(iorb) = max(0,inewocc(iorb)-1)
      end do

      ! Spin string loop:
      do indx=1,nstring

        ! Alpha index in full string space ...
        do i=1,nalfsing
          iaocc = iaccm(astr(i,indx))
          inewocc(iaocc) = inewocc(iaocc)+1
        end do
        iaind = indget_cvb(inewocc,nalf,norb,xalf)
        do i=1,nalfsing
          iaocc = iaccm(astr(i,indx))
          inewocc(iaocc) = inewocc(iaocc)-1
        end do

        ! Beta index in full string space ...
        do i=1,nbetsing
          ibocc = iaccm(bstr(i,indx))
          inewocc(ibocc) = inewocc(ibocc)+1
        end do
        ibind = indget_cvb(inewocc,nbet,norb,xbet)
        do i=1,nbetsing
          ibocc = iaccm(bstr(i,indx))
          inewocc(ibocc) = inewocc(ibocc)-1
        end do

        incrdet = incrdet+1
        idetavb(incrdet) = iaind
        idetbvb(incrdet) = ibind
      end do
    end do
    call mma_deallocate(astr)
    call mma_deallocate(bstr)
  end if
  ioff_nconf = ioff_nconf+nconfion(ion)
end do
if (debug) then
  write(u6,*) ' idetavb='
  write(u6,'(10i6)') idetavb
  write(u6,*) ' idetbvb='
  write(u6,'(10i6)') idetbvb
end if

call mma_deallocate(xalf)
call mma_deallocate(xbet)
call mma_deallocate(inewocc)
call mma_deallocate(iaccm)

return

end subroutine vbgenabdet_cvb
