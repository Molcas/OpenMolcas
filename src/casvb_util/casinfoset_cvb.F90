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

implicit real*8(a-h,o-z)
! ... Files/Hamiltonian available ...
logical, external :: valid_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "inpmod_cvb.fh"
#include "WrkSpc.fh"
#include "casinfo_cvb.fh"
logical hadinput
logical debug
data debug/.false./

if (debug) then
  write(6,*) ' casinfoset :'
  write(6,*) ' ------------'
  write(6,*) ' iorocc_c  :',iorocc_c
  write(6,*) ' iorclos_c :',iorclos_c
  write(6,*) ' iorcore_c :',iorcore_c
  write(6,*) ' nstsym_c  :',nstsym_c
  write(6,*) ' weight_c  :',weight_c
  write(6,*) ' istnel_c  :',istnel_c
  write(6,*) ' istsy_c   :',istsy_c
  write(6,*) ' istms2_c  :',istms2_c
  write(6,*) ' nstats_c  :',nstats_c
  write(6,*) ' strtint_c :',strtint_c
  write(6,*) ' strtci_c  :',strtci_c
  write(6,*) ' strtmo_c  :',strtmo_c
  write(6,*) ' iorocc_d  :',iorocc_d
  write(6,*) ' iorclos_d :',iorclos_d
  write(6,*) ' iorcore_d :',iorcore_d
  write(6,*) ' nstsym_d  :',nstsym_d
  write(6,*) ' weight_d  :',weight_d
  write(6,*) ' istnel_d  :',istnel_d
  write(6,*) ' istsy_d   :',istsy_d
  write(6,*) ' istms2_d  :',istms2_d
  write(6,*) ' nstats_d  :',nstats_d
  write(6,*) ' strtint_d :',strtint_d
  write(6,*) ' strtci_d  :',strtci_d
  write(6,*) ' strtmo_d  :',strtmo_d
end if
hadinput = .false.
do i=1,mxstsy_ci
  if (iorocc_d(i) /= -1) hadinput = .true.
  if (iorclos_d(i) /= -1) hadinput = .true.
  if (iorcore_d(i) /= -1) hadinput = .true.
end do
if (hadinput) then
  do i=1,mxstsy_ci
    if (iorocc_d(i) == -1) iorocc_d(i) = 0
    if (iorclos_d(i) == -1) iorclos_d(i) = 0
    if (iorcore_d(i) == -1) iorcore_d(i) = 0
  end do
else
  call imove_cvb(iorcore_c,iorcore_d,mxstsy_ci)
  call imove_cvb(iorclos_c,iorclos_d,mxstsy_ci)
  call imove_cvb(iorocc_c,iorocc_d,mxstsy_ci)
end if

mcore_d = 0
!  Ensure no negative number of orbitals:
do i=1,mxirrep
  iorclos_d(i) = iorclos_d(i)+iorcore_d(i)
  iorocc_d(i) = iorocc_d(i)+iorclos_d(i)
  mcore_d = mcore_d+iorclos_d(i)
end do

if (nstsym_d == 0) then
  nstsym_d = nstsym_c
  call imove_cvb(nstats_c,nstats_d,mxstsy_ci)
  call imove_cvb(istnel_c,istnel_d,mxstsy_ci)
  call imove_cvb(istsy_c,istsy_d,mxstsy_ci)
  call imove_cvb(istms2_c,istms2_d,mxstsy_ci)
  call fmove_cvb(weight_c,weight_d,mxstt_ci*mxstsy_ci)
  if (mcore_d /= mcore_c) then
!  Different number of core orbitals input -> assume NELTOT the same:
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

!  Set active space information
sum = zero
do i=1,nstsym_d
  do j=1,nstats_d(i)
    if (weight_d(j,i) < zero) then
      write(6,'(a,f10.4,i3,a,i1)') ' Fatal error: WEIGHT factor negative :',weight_d(j,i),j,'.',i
      call abend_cvb()
    end if
    sum = sum+weight_d(j,i)
  end do
end do
sum = one/sum
call dscal_(mxstt_ci*mxstsy_ci,sum,weight_d,1)
nel_d = -1
i2s_d = -1
call izero(isymv,mxirrep)
do i=1,nstsym_d
  do j=1,nstats_d(i)
    if (weight_d(j,i) > 1.d-20) then
      if ((nel_d /= -1) .and. (nel_d /= istnel_d(i))) then
        write(6,*) ' Fatal error: ELEC varies in WF cards!'
        call abend_cvb()
      end if
      if ((i2s_d /= -1) .and. (i2s_d /= istms2_d(i))) then
        write(6,*) ' Fatal error: SPIN varies in WF cards!'
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

call izero(ityp,mxorb_cvb)
mcore = 0
norb = 0
incr = 0
do i=1,mxirrep
  ioc = iorocc_d(i)-iorclos_d(i)
  norb = norb+ioc
  mcore = mcore+iorclos_d(i)
  do j=1,ioc
    ityp(j+incr) = i
  end do
  incr = incr+ioc
end do
!  Set NIRREP:
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
  write(6,*) ' Impossible numbers: active electrons :',nel
  write(6,*) '                     active orbitals  :',norb
  write(6,*) '                     total spin       :',dble(nalf-nbet)/two
  call abend_cvb()
end if
if (isym == 0) then
  write(6,*) ' WARNING: State symmetry not found - assuming A1.'
  isym = 1
  nsym = 1
  call izero(isymv,mxirrep)
  isymv(1) = 1
end if

return

end subroutine casinfoset_cvb
