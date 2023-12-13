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

subroutine casinfoset_cvb()

use casvb_global, only: iorclos_c, iorclos_d, iorcore_c, iorcore_d, iorocc_c, iorocc_d, istms2_c, istms2_d, istnel_c, istnel_d, &
                        istsy_c, istsy_d, isym, isymv, ityp, mcore_c, mcore_d, mxirrep, mxstsy_ci, nalf, nbet, nel, nirrep, noe, &
                        norb, nstats_c, nstats_d, nstsym_c, nstsym_d, nsym, strtci, strtint, strtmo, weight_c, weight_d
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, i2s_d, incr, ioc, irrep, is, isym_d, j, mcore, nel_d
real(kind=wp) :: rsum, strtci_c, strtci_d, strtint_c, strtint_d, strtmo_c, strtmo_d
logical(kind=iwp) :: hadinput
logical(kind=iwp), parameter :: debug = .false.
logical(kind=iwp), external :: valid_cvb ! ... Files/Hamiltonian available ...

! These were in a common block, but never initialized or used elsewhere
strtint_c = 0
strtci_c = 0
strtmo_c = 0
isym_d = 0

if (debug) then
  write(u6,*) ' casinfoset :'
  write(u6,*) ' ------------'
  write(u6,*) ' iorocc_c  :',iorocc_c
  write(u6,*) ' iorclos_c :',iorclos_c
  write(u6,*) ' iorcore_c :',iorcore_c
  write(u6,*) ' nstsym_c  :',nstsym_c
  write(u6,*) ' weight_c  :',weight_c
  write(u6,*) ' istnel_c  :',istnel_c
  write(u6,*) ' istsy_c   :',istsy_c
  write(u6,*) ' istms2_c  :',istms2_c
  write(u6,*) ' nstats_c  :',nstats_c
  write(u6,*) ' strtint_c :',strtint_c
  write(u6,*) ' strtci_c  :',strtci_c
  write(u6,*) ' strtmo_c  :',strtmo_c
  write(u6,*) ' iorocc_d  :',iorocc_d
  write(u6,*) ' iorclos_d :',iorclos_d
  write(u6,*) ' iorcore_d :',iorcore_d
  write(u6,*) ' nstsym_d  :',nstsym_d
  write(u6,*) ' weight_d  :',weight_d
  write(u6,*) ' istnel_d  :',istnel_d
  write(u6,*) ' istsy_d   :',istsy_d
  write(u6,*) ' istms2_d  :',istms2_d
  write(u6,*) ' nstats_d  :',nstats_d
  write(u6,*) ' strtint_d :',strtint_d
  write(u6,*) ' strtci_d  :',strtci_d
  write(u6,*) ' strtmo_d  :',strtmo_d
end if
hadinput = .false.
do i=1,mxstsy_ci
  if ((iorocc_d(i) /= -1) .or. (iorclos_d(i) /= -1) .or. (iorcore_d(i) /= -1)) then
    hadinput = .true.
    exit
  end if
end do
if (hadinput) then
  do i=1,mxstsy_ci
    if (iorocc_d(i) == -1) iorocc_d(i) = 0
    if (iorclos_d(i) == -1) iorclos_d(i) = 0
    if (iorcore_d(i) == -1) iorcore_d(i) = 0
  end do
else
  iorcore_d(:) = iorcore_c(:)
  iorclos_d(:) = iorclos_c(:)
  iorocc_d(:) = iorocc_c(:)
end if

mcore_d = 0
! Ensure no negative number of orbitals:
iorclos_d(:) = iorclos_d(:)+iorcore_d(:)
iorocc_d(:) = iorocc_d(:)+iorclos_d(:)
mcore_d = sum(iorclos_d(:))

if (nstsym_d == 0) then
  nstsym_d = nstsym_c
  nstats_d(:) = nstats_c(:)
  istnel_d(:) = istnel_c(:)
  istsy_d(:) = istsy_c(:)
  istms2_d(:) = istms2_c(:)
  weight_d(:,:) = weight_c(:,:)
  if (mcore_d /= mcore_c) then
    ! Different number of core orbitals input -> assume NELTOT the same:
    do i=1,mxstsy_ci
      if (istnel_d(i) /= 0) istnel_d(i) = istnel_d(i)+2*(mcore_c-mcore_d)
    end do
  end if
end if
strtint_d = strtint
strtmo_d = strtmo
strtci_d = strtci
if (.not. valid_cvb(strtint_d)) strtint_d = strtint_c
if (.not. valid_cvb(strtmo_d)) strtmo_d = strtmo_c
if (.not. valid_cvb(strtci_d)) strtci_d = strtci_c
strtint = strtint_d
strtmo = strtmo_d
strtci = strtci_d

! Set active space information
rsum = Zero
do i=1,nstsym_d
  do j=1,nstats_d(i)
    if (weight_d(j,i) < Zero) then
      write(u6,'(a,f10.4,i3,a,i1)') ' Fatal error: WEIGHT factor negative :',weight_d(j,i),j,'.',i
      call abend_cvb()
    end if
    rsum = rsum+weight_d(j,i)
  end do
end do
weight_d(:,:) = weight_d(:,:)/rsum
nel_d = -1
i2s_d = -1
isymv(:) = 0
do i=1,nstsym_d
  do j=1,nstats_d(i)
    if (weight_d(j,i) > 1.0e-20_wp) then
      if ((nel_d /= -1) .and. (nel_d /= istnel_d(i))) then
        write(u6,*) ' Fatal error: ELEC varies in WF cards!'
        call abend_cvb()
      end if
      if ((i2s_d /= -1) .and. (i2s_d /= istms2_d(i))) then
        write(u6,*) ' Fatal error: SPIN varies in WF cards!'
        call abend_cvb()
      end if
      nel_d = istnel_d(i)
      i2s_d = istms2_d(i)
      isym_d = istsy_d(i)
      isymv(isym_d) = 1
      exit
    end if
  end do
end do
nsym = 0
do is=1,mxirrep
  if (isymv(is) == 1) nsym = nsym+1
end do

nel = nel_d
isym = isym_d

ityp(:) = 0
mcore = 0
norb = 0
incr = 0
do i=1,mxirrep
  ioc = iorocc_d(i)-iorclos_d(i)
  norb = norb+ioc
  mcore = mcore+iorclos_d(i)
  ityp(incr+1:incr+ioc) = i
  incr = incr+ioc
end do
! Set NIRREP:
nirrep = 1
do irrep=1,mxirrep
  if ((iorcore_d(irrep) > 0) .or. (iorclos_d(irrep) > 0) .or. (iorocc_d(irrep) > 0)) nirrep = irrep
end do
if (nirrep == 3) nirrep = 4
if (nirrep > 4) nirrep = 8

noe = max(norb,nel)
nbet = (nel-i2s_d)/2
nalf = nel-nbet
! Basic checks
if ((nel < 0) .or. (norb < 0) .or. (i2s_d < 0) .or. (nel > 2*norb) .or. (mod(nel,2) /= mod(i2s_d,2))) then
  write(u6,*) ' Impossible numbers: active electrons :',nel
  write(u6,*) '                     active orbitals  :',norb
  write(u6,*) '                     total spin       :',real(nalf-nbet,kind=wp)*Half
  call abend_cvb()
end if
if (isym == 0) then
  write(u6,*) ' WARNING: State symmetry not found - assuming A1.'
  isym = 1
  nsym = 1
  isymv(1) = 1
  isymv(2:) = 0
end if

return

end subroutine casinfoset_cvb
