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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_Drv(irc,EMP2,CMO,EOcc,EVir)
! Francesco Aquilante, May 2007.
!
! Purpose: driver for computing the Scaled Opposite-Spin (SOS)
!          MP2 energy correction EMP2
!          using Cholesky (or RI) representation for the
!          two-electron integrals.
!          Input must have been processed and MO coefficients
!          and orbital energies must be passed as arguments.
!
! Notes:
!
!   - all MO Cholesky vector files generated here are deleted before
!     exit, except for error terminations (i.e. no cleanup actions
!     are taken!)

use Cholesky, only: LuPri, nSym, NumCho
use ChoMP2, only: nMP2Vec, nT1am, set_cd_thr, ThrMP2, Verbose
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Five
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(in) :: CMO(*), EOcc(*), EVir(*)
integer(kind=iwp) :: iSym, lDiag, nSym_Sav
real(kind=wp) :: CPUDec1, CPUDec2, CPUEnr1, CPUEnr2, CPUIni1, CPUIni2, CPUTot1, CPUTot2, CPUTra1, CPUTra2, Diff, Dum, FracMem, &
                 WallDec1, WallDec2, WallEnr1, WallEnr2, WallIni1, WallIni2, WallTot1, WallTot2, WallTra1, WallTra2
logical(kind=iwp) :: Delete
real(kind=wp), allocatable :: Diag(:)
integer(kind=iwp), parameter :: iFmt = 0
real(kind=wp), parameter :: Chk_Mem_ChoMP2 = 0.123456789_wp, Tol = 1.0e-15_wp
logical(kind=iwp), parameter :: Delete_def = .true.
character(len=*), parameter :: SecNam = 'Cho_SOSmp2_Drv'
real(kind=wp), external :: ddot_

#ifdef _DEBUGPRINT_
Verbose = .true.
#endif
if (Verbose) then
  call CWTime(CPUTot1,WallTot1)
end if

! Initializations.
! ----------------

irc = 0

EMP2 = Zero

if (Verbose) then
  call CWTime(CPUIni1,WallIni1)
end if

Dum = Chk_Mem_ChoMP2

FracMem = Zero ! no buffer allocated
call Cho_X_Init(irc,FracMem)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_X_Init returned ',irc
  call SysAbendMsg(SecNam,'Cholesky initialization error',' ')
end if

call Cho_SOSmp2_Setup(irc)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_SOSmp2_Setup returned ',irc
  call finalize()
  return
end if

if (Verbose) then
  call Cho_SOSmp2_Setup_Prt(irc)
  if (irc /= 0) then
    write(u6,*) SecNam,': Cho_SOSmp2_Setup_Prt returned ',irc
    call finalize()
    return
  end if
  call CWTime(CPUIni2,WallIni2)
  call Cho_PrtTim('Cholesky SOS-MP2 initialization',CPUIni2,CPUIni1,WallIni2,WallIni1,iFmt)
end if

! Transform Cholesky vectors directly from reduced set to MO
! representation. Result vectors are stored on disk.
! Compute also the (ai|ai)^2 diagonal here.
! ----------------------------------------------------------

if (Verbose) then
  call CWTime(CPUTra1,WallTra1)
end if
lDiag = nT1am(1)
do iSym=2,nSym
  lDiag = lDiag+nT1am(iSym)
end do
call mma_allocate(Diag,lDiag,label='Diag')

call ChoMP2_TraDrv(irc,CMO,Diag,.true.)
if (irc /= 0) then
  write(u6,*) SecNam,': ChoMP2_TraDrv returned ',irc
  call finalize()
  return
end if

! Squaring each diagonal element
! ------------------------------
Diag(:) = Diag(:)**2
if (set_cd_thr) ThrMP2 = ddot_(lDiag,[One],0,Diag,1)/(Five*lDiag)

if (Verbose) then
  call CWTime(CPUTra2,WallTra2)
  call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,CPUTra1,WallTra2,WallTra1,iFmt)
end if

! Finalize Cholesky info (to release memory).
! Retain essential info: LuPri, nSym, and NumCho(*).
! --------------------------------------------------

nSym_Sav = nSym
nMP2Vec(1:nSym) = NumCho(1:nSym)

call Cho_X_Final(irc)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_X_Final returned ',irc
  call finalize()
  return
end if

LuPri = u6
nSym = nSym_Sav
NumCho(1:nSym) = nMP2Vec(1:nSym)

! Decompose M(ai,bj) = (ai|bj)^2 .
! Set number of vectors to be used in energy calculation.
! -------------------------------------------------------

if (Verbose) then
  call CWTime(CPUDec1,WallDec1)
end if
Delete = Delete_def ! delete transf. vector files after dec.
call Cho_SOSmp2_DecDrv(irc,Delete,Diag)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_SOSmp2_DecDrv returned ',irc
  call SysAbendMsg(SecNam,'SOS-MP2 decomposition failed!',' ')
end if
if (Verbose) then
  call CWTime(CPUDec2,WallDec2)
  call Cho_PrtTim('Cholesky SOS-MP2 decomposition',CPUDec2,CPUDec1,WallDec2,WallDec1,iFmt)
end if
call mma_deallocate(Diag)

! Compute SOS-MP2 energy correction.
! ----------------------------------

if (Verbose) then
  call CWTime(CPUEnr1,WallEnr1)
end if
Delete = Delete_def
call Cho_SOSmp2_Energy(irc,EMP2,EOcc,EVir,Delete)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_SOSmp2_Energy returned ',irc
  call finalize
  return
end if
if (Verbose) then
  call CWTime(CPUEnr2,WallEnr2)
  call Cho_PrtTim('Cholesky SOS-MP2 energy',CPUEnr2,CPUEnr1,WallEnr2,WallEnr1,iFmt)
end if

! Exit.
! -----

call finalize()

return

contains

subroutine finalize()
  Diff = abs(Dum-Chk_Mem_ChoMP2)
  if (Diff > Tol) then
    write(u6,*) SecNam,': Memory Boundary Error!'
    if (irc == 0) irc = -9999
  end if
  if (Verbose) then
    call CWTime(CPUTot2,WallTot2)
    call Cho_PrtTim('Cholesky SOS-MP2',CPUTot2,CPUTot1,WallTot2,WallTot1,iFmt)
  end if
end subroutine finalize

end subroutine Cho_SOSmp2_Drv
