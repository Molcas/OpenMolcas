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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ChoLSOSMP2_Energy(irc,EMP2,EOcc,EVir,Sorted,DelOrig)
!
! Thomas Bondo Pedersen, December 2012.
!
! Compute Laplace-SOS-MP2 energy.

use stdalloc

implicit none
integer irc
real*8 EMP2
real*8 EOcc(*)
real*8 EVir(*)
logical Sorted
logical DelOrig

#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"
character(len=17), parameter :: SecNam = 'ChoLSOSMP2_Energy'
integer CheckDenomRange
integer, external :: TestMinimaxLaplace
#ifdef _DEBUGPRINT_
logical, parameter :: Debug = .true.
#else
logical, parameter :: Debug = .false.
#endif
logical Verb, FermiShift
integer l_w
real*8 ELOMO, EHOMO
real*8 ELUMO, EHUMO
real*8 EFermi
real*8 xmin, xmax
integer iSym
integer i
real*8, allocatable :: W(:), T(:)

!================
! Initializations
!================

call Untested('Laplace-SOS-MP2')

! init return code
irc = 0

! init flag for Fermi shift done
FermiShift = .false.

! check that Laplace is requested
if (.not. Laplace) then
  call WarningMessage(1,SecNam//' was called - but this is not a Laplace calculation!')
  call xFlush(6)
  irc = -1
  return
end if

! Debug: test minimax Laplace grid generation
if (Debug) then
  Verb = .false.
  irc = TestMinimaxLaplace(1.0d-7,Verb)
  if (irc /= 0) then
    call WarningMessage(2,SecNam//': error detected in numerical Laplace transformation')
    write(6,'(A,I6)') 'irc=',irc
    call Abend()
  end if
end if

!================================================
! Parameters for numerical Laplace transformation
!================================================

! get max and min orbital energies
ELOMO = 0.0d0
EHOMO = 0.0d0
ELUMO = 0.0d0
EHUMO = 0.0d0
i = 0
do iSym=1,nSym
  if (nOcc(iSym) > 0) then
    if (i == 0) then
      i = 1
      ELOMO = EOcc(iOcc(iSym)+1)
      EHOMO = EOcc(iOcc(iSym)+nOcc(iSym))
    else
      ELOMO = min(ELOMO,EOcc(iOcc(iSym)+1))
      EHOMO = max(EHOMO,EOcc(iOcc(iSym)+nOcc(iSym)))
    end if
  end if
end do
if (i == 0) then
  call WarningMessage(2,SecNam//': unable to determine LOMO,HOMO')
  call Abend()
end if
i = 0
do iSym=1,nSym
  if (nVir(iSym) > 0) then
    if (i == 0) then
      i = 1
      ELUMO = EVir(iVir(iSym)+1)
      EHUMO = EVir(iVir(iSym)+nVir(iSym))
    else
      ELUMO = min(ELUMO,EVir(iVir(iSym)+1))
      EHUMO = max(EHUMO,EVir(iVir(iSym)+nVir(iSym)))
    end if
  end if
end do
if (i == 0) then
  call WarningMessage(2,SecNam//': unable to determine LUMO,HUMO')
  call Abend()
end if
!-tbp:
write(6,*) 'ELOMO,EHOMO=',ELOMO,EHOMO
write(6,*) 'ELUMO,EHUMO=',ELUMO,EHUMO

! compute "Fermi energy" as the midpoint between HOMO and LUMO.
EFermi = 0.5d0*(EHOMO+ELUMO)

! translate orbital energy origin to EFermi
do iSym=1,nSym
  do i=1,nOcc(iSym)
    EOcc(iOcc(iSym)+i) = EOcc(iOcc(iSym)+i)-EFermi
  end do
  do i=1,nVir(iSym)
    EVir(iVir(iSym)+i) = EVir(iVir(iSym)+i)-EFermi
  end do
end do
ELOMO = ELOMO-EFermi
EHOMO = EHOMO-EFermi
ELUMO = ELUMO-EFermi
EHUMO = EHUMO-EFermi
FermiShift = .true.

! compute range of orbital energy denominator
xmin = 2.0d0*(ELUMO-EHOMO)
xmax = 2.0d0*(EHUMO-ELOMO)
! Debug: check range
if (Debug) then
  irc = CheckDenomRange(xmin,xmax,nSym,EOcc,Evir,iOcc,nOcc,iVir,nVir)
  if (irc /= 0) then
    call WarningMessage(2,SecNam//': error detected in orbital energy denominator range')
    write(6,'(A,I6)') 'irc=',irc
    call Abend()
  end if
end if

! get weights and grid points for numerical Laplace transform
if (Laplace_nGridPoints == 0) then
  l_w = Laplace_mGridPoints
else
  l_w = Laplace_nGridPoints
end if
call mma_allocate(W,l_w,Label='W')
call mma_allocate(T,l_w,Label='T')
call MinimaxLaplace(Verbose,Laplace_nGridPoints,xmin,xmax,l_w,W,T,irc)
if (irc /= 0) then
  write(6,'(A,A,I6)') SecNam,': MinimaxLaplace returned',irc
  irc = 1
  Go To 1 ! exit after cleanup actions
end if

!==================================
! Compute SOS-MP2 energy correction
!==================================

if (Sorted) then
  call ChoLSOSMP2_Energy_Srt(Laplace_nGridPoints,W,T,EOcc,EVir,DelOrig,EMP2,irc)
  if (irc /= 0) then
    write(6,'(A,A,I6)') SecNam,': ChoLSOSMP2_Energy_Srt returned',irc
    Go To 1 ! exit
  end if
else
  if (nBatch == 1) then
    call ChoLSOSMP2_Energy_Fll(Laplace_nGridPoints,W,T,EOcc,EVir,DelOrig,EMP2,irc)
    if (irc /= 0) then
      write(6,'(A,A,I6)') SecNam,': ChoLSOSMP2_Energy_Fll returned',irc
      Go To 1 ! exit after cleanup
    end if
  else
    call WarningMessage(1,SecNam//': unsorted case not implemented')
    irc = -2
    Go To 1 ! exit after cleanup
  end if
end if

!========
! Cleanup
!========
1 continue ! errors jump to this point

! translate orbital energy origin from Fermi back to original
if (FermiShift) then
  do iSym=1,nSym
    do i=1,nOcc(iSym)
      EOcc(iOcc(iSym)+i) = EOcc(iOcc(iSym)+i)+EFermi
    end do
    do i=1,nVir(iSym)
      EVir(iVir(iSym)+i) = EVir(iVir(iSym)+i)+EFermi
    end do
  end do
end if

! deallocations
call mma_deallocate(T)
call mma_deallocate(W)

end subroutine ChoLSOSMP2_Energy
