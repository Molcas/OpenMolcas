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

use Cholesky, only: nSym
use ChoMP2, only: iOcc, iVir, Laplace, Laplace_mGridPoints, Laplace_nGridPoints, nBatch, nOcc, nVir, Verbose
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(inout) :: EOcc(*), EVir(*)
logical(kind=iwp), intent(in) :: Sorted, DelOrig
integer(kind=iwp) :: CheckDenomRange, i, iSym, l_w
logical(kind=iwp) :: FermiShift, Verb
real(kind=wp) :: EFermi, EHOMO, EHUMO, ELOMO, ELUMO, xmax, xmin
real(kind=wp), allocatable :: W(:), T(:)
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: Debug = _DBG_
character(len=*), parameter :: SecNam = 'ChoLSOSMP2_Energy'
integer(kind=iwp), external :: TestMinimaxLaplace

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
  call xFlush(u6)
  irc = -1
  return
end if

! Debug: test minimax Laplace grid generation
if (Debug) then
  Verb = .false.
  irc = TestMinimaxLaplace(1.0e-7_wp,Verb)
  if (irc /= 0) then
    call WarningMessage(2,SecNam//': error detected in numerical Laplace transformation')
    write(u6,'(A,I6)') 'irc=',irc
    call Abend()
  end if
end if

!================================================
! Parameters for numerical Laplace transformation
!================================================

! get max and min orbital energies
ELOMO = Zero
EHOMO = Zero
ELUMO = Zero
EHUMO = Zero
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
write(u6,*) 'ELOMO,EHOMO=',ELOMO,EHOMO
write(u6,*) 'ELUMO,EHUMO=',ELUMO,EHUMO

! compute "Fermi energy" as the midpoint between HOMO and LUMO.
EFermi = Half*(EHOMO+ELUMO)

! translate orbital energy origin to EFermi
do iSym=1,nSym
  EOcc(iOcc(iSym)+1:iOcc(iSym)+nOcc(iSym)) = EOcc(iOcc(iSym)+1:iOcc(iSym)+nOcc(iSym))-EFermi
  EVir(iVir(iSym)+1:iVir(iSym)+nVir(iSym)) = EVir(iVir(iSym)+1:iVir(iSym)+nVir(iSym))-EFermi
end do
ELOMO = ELOMO-EFermi
EHOMO = EHOMO-EFermi
ELUMO = ELUMO-EFermi
EHUMO = EHUMO-EFermi
FermiShift = .true.

! compute range of orbital energy denominator
xmin = Two*(ELUMO-EHOMO)
xmax = Two*(EHUMO-ELOMO)
! Debug: check range
if (Debug) then
  irc = CheckDenomRange(xmin,xmax,nSym,EOcc,Evir,iOcc,nOcc,iVir,nVir)
  if (irc /= 0) then
    call WarningMessage(2,SecNam//': error detected in orbital energy denominator range')
    write(u6,'(A,I6)') 'irc=',irc
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
  write(u6,'(A,A,I6)') SecNam,': MinimaxLaplace returned',irc
  irc = 1
  ! exit after cleanup actions
  call Finish_this()
  return
end if

!==================================
! Compute SOS-MP2 energy correction
!==================================

if (Sorted) then
  call ChoLSOSMP2_Energy_Srt(Laplace_nGridPoints,W,T,EOcc,EVir,DelOrig,EMP2,irc)
  if (irc /= 0) then
    write(u6,'(A,A,I6)') SecNam,': ChoLSOSMP2_Energy_Srt returned',irc
    call Finish_this()
    return
  end if
else
  if (nBatch == 1) then
    call ChoLSOSMP2_Energy_Fll(Laplace_nGridPoints,W,T,EOcc,EVir,DelOrig,EMP2,irc)
    if (irc /= 0) then
      write(u6,'(A,A,I6)') SecNam,': ChoLSOSMP2_Energy_Fll returned',irc
      ! exit after cleanup
      call Finish_this()
      return
    end if
  else
    call WarningMessage(1,SecNam//': unsorted case not implemented')
    irc = -2
    ! exit after cleanup
    call Finish_this()
    return
  end if
end if

call Finish_this()

contains

!========
! Cleanup
!========
subroutine Finish_this()

  integer(kind=iwp) :: iSym

  ! translate orbital energy origin from Fermi back to original
  if (FermiShift) then
    do iSym=1,nSym
      EOcc(iOcc(iSym)+1:iOcc(iSym)+nOcc(iSym)) = EOcc(iOcc(iSym)+1:iOcc(iSym)+nOcc(iSym))+EFermi
      EVir(iVir(iSym)+1:iVir(iSym)+nVir(iSym)) = EVir(iVir(iSym)+1:iVir(iSym)+nVir(iSym))+EFermi
    end do
  end if

  ! deallocations
  call mma_deallocate(T)
  call mma_deallocate(W)

end subroutine Finish_this

end subroutine ChoLSOSMP2_Energy
