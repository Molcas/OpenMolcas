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

subroutine cnfcheck2_cvb(iconfs,nconf1,nel1)

use casvb_global, only: noe, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: nconf1, iconfs(noe,nconf1)
integer(kind=iwp), intent(in) :: nel1
integer(kind=iwp) :: i, iconf, ii, iorb, jconf, nsum
logical(kind=iwp) :: found, locc, locc_only, lorbs, lorbs_only
integer(kind=iwp), allocatable :: iocc(:)

if (nconf1 == 0) then
  ! Special case -- iconfs will be ok and adhere to occ no definition.
  nconf1 = 1
  return
end if

! Perform basic checks of configurations:
! First determine if a consistent definition (orb list or
! occ numbers) has been used for the configurations.
!
! We do this in two parses. locc_only/lorbs_only are set in
! first run through if *any* iconf unamiguously adheres to
! one of the definitions.
!
! Second run through if necessary we check again in case iconfs
! contains *both* types of definitions.

call mma_allocate(iocc,noe,label='iocc')

locc_only = .false.
lorbs_only = .false.
do iconf=1,nconf1
  ! Consistency with occ no definition?
  locc = .true.
  do i=norb+1,noe
    if (iconfs(i,iconf) /= 0) then
      locc = .false.
      exit
    end if
  end do
  nsum = 0
  do i=1,norb
    if ((iconfs(i,iconf) < 0) .or. (iconfs(i,iconf) > 2)) locc = .false.
    nsum = nsum+iconfs(i,iconf)
  end do
  if (nsum /= nel1) locc = .false.
  ! Consistency with orb list definition?
  lorbs = .true.
  do i=nel1+1,noe
    if (iconfs(i,iconf) /= 0) then
      lorbs = .false.
      exit
    end if
  end do
  iocc(1:norb) = 0
  do i=1,nel1
    if ((iconfs(i,iconf) >= 1) .and. (iconfs(i,iconf) <= norb)) then
      iocc(iconfs(i,iconf)) = iocc(iconfs(i,iconf))+1
    else
      lorbs = .false.
    end if
  end do
  do i=1,norb
    if (iocc(i) > 2) then
      lorbs = .false.
      exit
    end if
  end do

  if (locc .and. (.not. lorbs)) then
    locc_only = .true.
  else if (lorbs .and. (.not. locc)) then
    lorbs_only = .true.
  else if ((.not. lorbs) .and. (.not. locc)) then
    write(u6,*) ' Illegal configuration read ',iconf
    write(u6,*) iconfs(1:noe,iconf)
    call abend_cvb()
  end if
end do

locc = locc_only
lorbs = lorbs_only
do iconf=1,nconf1
  if (locc_only .and. lorbs_only) then
    ! Check again ...
    ! Consistency with occ no definition?
    locc = .true.
    do i=norb+1,noe
      if (iconfs(i,iconf) /= 0) then
        locc = .false.
        exit
      end if
    end do
    nsum = 0
    do i=1,norb
      if ((iconfs(i,iconf) < 0) .or. (iconfs(i,iconf) > 2)) locc = .false.
      nsum = nsum+iconfs(i,iconf)
    end do
    if (nsum /= nel1) locc = .false.
    ! Consistency with orb list definition?
    lorbs = .true.
    do i=nel1+1,noe
      if (iconfs(i,iconf) /= 0) then
        lorbs = .false.
        exit
      end if
    end do
    iocc(1:norb) = 0
    do i=1,nel1
      if ((iconfs(i,iconf) >= 1) .and. (iconfs(i,iconf) <= norb)) then
        iocc(iconfs(i,iconf)) = iocc(iconfs(i,iconf))+1
      else
        lorbs = .false.
      end if
    end do
    do i=1,norb
      if (iocc(i) > 2) then
        lorbs = .false.
        exit
      end if
    end do
  end if
  if (locc .and. lorbs) then
    ! Comment out following 5 lines if default should be occ no definition:
    iocc(1:norb) = 0
    do i=1,nel1
      iocc(iconfs(i,iconf)) = iocc(iconfs(i,iconf))+1
    end do
    iconfs(1:norb,iconf) = iocc(1:norb)
    if (noe-norb > 0) iconfs(norb+1:,iconf) = 0
  else if (lorbs) then
    iocc(1:norb) = 0
    do i=1,nel1
      iocc(iconfs(i,iconf)) = iocc(iconfs(i,iconf))+1
    end do
    iconfs(1:norb,iconf) = iocc(1:norb)
    if (noe-norb > 0) iconfs(norb+1:,iconf) = 0
  end if
  if (iconf <= 500) then
    ! Test for repeated configurations:
    do jconf=1,iconf-1
      found = .false.
      do iorb=1,norb
        if (iconfs(iorb,iconf) /= iconfs(iorb,jconf)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        write(u6,'(/,a,2i4)') ' Fatal error - spatial VB configuration repeated :',jconf,iconf
        write(u6,'(i8,a,20i3)') jconf,'   =>  ',(iocc(ii),ii=1,norb)
        write(u6,'(i8,a,20i3)') iconf,'   =>  ',(iocc(ii),ii=1,norb)
        call abend_cvb()
      end if
    end do
  end if
end do

call mma_deallocate(iocc)

return

end subroutine cnfcheck2_cvb
