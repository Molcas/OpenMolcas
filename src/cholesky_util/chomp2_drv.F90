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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Drv(irc,EMP2,CMO,EOcc,EVir,Dab,Dii)
!
! Thomas Bondo Pedersen, October 2004.
!
! Purpose: driver for computing the MP2 energy correction EMP2
!          using Cholesky decomposed two-electron integrals.
!          Input must have been processed and MO coefficients
!          and orbital energies must be passed as arguments.
!
! Notes:
!
!   - all MO Cholesky vector files generated here are deleted before
!     exit, except for error terminations (i.e. no cleanup actions
!     are taken!)

use Cholesky, only: LuPri, nSym, NumCho
use ChoMP2, only: DecoMP2, DoDens, DoFNO, DoGrdt, EFrozT, EMP2_dens, EOccuT, EVirtT, l_Dii, Laplace, nBatch, nMoMo, nMP2Vec, &
                  nT1am, SOS_mp2, Verbose
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(inout) :: EOcc(*), EVir(*), Dab(*), Dii(*)
integer(kind=iwp) :: lDiag, nSym_Sav
real(kind=wp) :: CPUDab1, CPUDab2, CPUDec1, CPUDec2, CPUEnr1, CPUEnr2, CPUIni1, CPUIni2, CPUSrt1, CPUSrt2, CPUTot1, CPUTot2, &
                 CPUTra1, CPUTra2, FracMem, WallDab1, WallDab2, WallDec1, WallDec2, WallEnr1, WallEnr2, WallIni1, WallIni2, &
                 WallSrt1, WallSrt2, WallTot1, WallTot2, WallTra1, WallTra2
logical(kind=iwp) :: Delete, DoSort
real(kind=wp), allocatable :: Check(:), Diag(:)
integer(kind=iwp), parameter :: iFmt = 0
real(kind=wp), parameter :: Chk_Mem_ChoMP2 = 0.123456789_wp, Tol = 1.0e-15_wp
logical(kind=iwp), parameter :: Delete_def = .true.
character(len=*), parameter :: SecNam = 'ChoMP2_Drv'

#ifdef _DEBUGPRINT_
Verbose = .true.
#endif
if (Verbose) call CWTime(CPUTot1,WallTot1)

! Initializations.
! ----------------

irc = 0

EMP2 = Zero
if (DoDens) EMP2_dens = Zero

if (Verbose) call CWTime(CPUIni1,WallIni1)

call mma_allocate(Check,1,Label='Check')
Check(1) = Chk_Mem_ChoMP2

FracMem = Zero ! no buffer allocated
call Cho_X_Init(irc,FracMem)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_X_Init returned ',irc
  call SysAbendMsg(SecNam,'Cholesky initialization error',' ')
end if

call ChoMP2_Setup(irc)
if (irc /= 0) then
  write(u6,*) SecNam,': ChoMP2_Setup returned ',irc
  call Finish_this()
  return
end if

if (DoDens) then
  !write(u6,*) 'Run ChoMP2g_setup'
  call ChoMP2g_Setup(irc,EOcc,EVir)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2g_Setup returned ',irc
    call Finish_this()
    return
  end if
end if

if (Verbose) then
  call ChoMP2_Setup_Prt(irc)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_Setup_Prt returned ',irc
    call Finish_this()
    return
  end if
  call CWTime(CPUIni2,WallIni2)
  call Cho_PrtTim('Cholesky MP2 initialization',CPUIni2,CPUIni1,WallIni2,WallIni1,iFmt)
end if

! Transform Cholesky vectors directly from reduced set to MO
! representation. Result vectors are stored on disk.
! If decomposition of (ai|bj) is requested, compute also the
! (ai|ai) diagonal here.
! ----------------------------------------------------------

if (Verbose) call CWTime(CPUTra1,WallTra1)
if (DecoMP2) then
  lDiag = sum(nT1am(1:nSym))
else if (DoDens) then
  lDiag = sum(nMoMo(1:nSym,6))
else
  lDiag = 1
end if
call mma_allocate(Diag,lDiag,Label='Diag')

if (.not. DoDens) then
  call ChoMP2_TraDrv(irc,CMO,Diag,DecoMP2)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_TraDrv returned ',irc
    call Finish_this()
    return
  end if
  if (Verbose) then
    call CWTime(CPUTra2,WallTra2)
    call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,CPUTra1,WallTra2,WallTra1,iFmt)
  end if
else if (DoDens) then
  call ChoMP2g_TraDrv(irc,CMO,Diag,DecoMP2)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2g_TraDrv returned ',irc
    call Finish_this()
    return
  end if
  if (Verbose) then
    call CWTime(CPUTra2,WallTra2)
    call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,CPUTra1,WallTra2,WallTra1,iFmt)
  end if
end if

! Finalize Cholesky info (to release memory).
! Retain essential info: LuPri, nSym, and NumCho(*).
! --------------------------------------------------

nSym_Sav = nSym
nMP2Vec(1:nSym) = NumCho(1:nSym)

call Cho_X_Final(irc)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_X_Final returned ',irc
  call Finish_this()
  return
end if

LuPri = u6
nSym = nSym_Sav
NumCho(1:nSym) = nMP2Vec(1:nSym)

! Decompose (ai|bj) integrals, if requested.
! Set number of vectors to be used in energy calculation.
! -------------------------------------------------------

if (DecoMP2) then
  if (Verbose) call CWTime(CPUDec1,WallDec1)
  Delete = Delete_def ! delete transf. vector files after dec.
  call ChoMP2_DecDrv(irc,Delete,Diag,'Integrals')
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_DecDrv returned ',irc
    call SysAbendMsg(SecNam,'MP2 decomposition failed!',' ')
  end if
  if (Verbose) then
    call CWTime(CPUDec2,WallDec2)
    call Cho_PrtTim('Cholesky MP2 decomposition',CPUDec2,CPUDec1,WallDec2,WallDec1,iFmt)
  end if
else if (DoDens) then
  if (Verbose) call CWTime(CPUDec1,WallDec1)
  call ChoMP2g_AmpDiag(irc,Diag,EOcc,EVir)
  Delete = .false.    ! do not delete transf. vectors.
  call ChoMP2_DecDrv(irc,Delete,Diag,'Amplitudes')
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_DecDrv returned ',irc
    call SysAbendMsg(SecNam,'MP2 decomposition failed!',' ')
  end if
  if (Verbose) then
    call CWTime(CPUDec2,WallDec2)
    call Cho_PrtTim('Cholesky MP2 decomposition',CPUDec2,CPUDec1,WallDec2,WallDec1,iFmt)
  end if
else
  nMP2Vec(1:nSym) = NumCho(1:nSym)
end if
call mma_deallocate(Diag)

! Presort Cholesky vectors if needed.
! -----------------------------------

DoSort = nBatch > 1
if (DoSort .and. (.not. DoDens)) then
  if (Verbose) call CWTime(CPUSrt1,WallSrt1)
  Delete = Delete_def
  call ChoMP2_SrtDrv(irc,Delete)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_SrtDrv returned ',irc
    if (Delete) then ! full vectors not available
      call SysAbendMsg(SecNam,'MP2 presort failed!',' ')
    else
      write(u6,*) SecNam,': trying to use full vectors instead...'
    end if
    DoSort = .false.
  end if
  if (Verbose) then
    call CWTime(CPUSrt2,WallSrt2)
    call Cho_PrtTim('Cholesky MP2 presort',CPUSrt2,CPUSrt1,WallSrt2,WallSrt1,iFmt)
  end if
end if

! FNO section: MP2 pseudodensity
! ------------------------------

if (DoFNO .and. (.not. DoDens)) then
  if (Verbose) call CWTime(CPUDab1,WallDab1)
  Delete = Delete_def
  call ChoMP2_FNO(irc,Dab,Dii,EOcc,EVir,DoSort,Delete)
  Dii(1:l_Dii) = -Dii(1:l_Dii)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_FNO returned ',irc
    call Finish_this()
    return
  end if
  if (Verbose) then
    call CWTime(CPUDab2,WallDab2)
    call Cho_PrtTim('Cholesky MP2 FNO section ',CPUDab2,CPUDab1,WallDab2,WallDab1,iFmt)
  end if
  call Finish_this()
  return
end if

! Compute MP2 Density and energy correction.
! ------------------------------------------
if (DoDens) then
  if (Verbose) call CWTime(CPUEnr1,WallEnr1)
  Delete = .false.
  call ChoMP2g_DensDrv(irc,EOccuT,EVirtT,EFrozT,CMO)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2g_DensDrv returned ',irc
    call Finish_this()
    return
  end if
end if

! Compute some matrices for Mp2-gradients
! ---------------------------------------
if (DoGrdt) then
  if (Verbose) call CWTime(CPUEnr1,WallEnr1)
  call ChoMP2g_GradSetup(irc,CMO)
  if (irc /= 0) then
    write(u6,*) SecNam,':ChoMP2g_GradSetup returned ',irc
    call Finish_this()
    return
  end if
  if (Verbose) then
    call CWTime(CPUEnr2,WallEnr2)
    call Cho_PrtTim('Cholesky Grad setup',CPUEnr2,CPUEnr1,WallEnr2,WallEnr1,iFmt)
  end if
end if

! Compute MP2 energy correction.
! ------------------------------

if (Verbose) call CWTime(CPUEnr1,WallEnr1)
Delete = Delete_def
if (Laplace .and. SOS_MP2) then
  call ChoLSOSMP2_Energy(irc,EMP2,EOcc,EVir,DoSort,Delete)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoLSOSMP2_Energy returned ',irc
    call Finish_this()
    return
  end if
else
  call ChoMP2_Energy(irc,EMP2,EOcc,EVir,DoSort,Delete)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_Energy returned ',irc
    call Finish_this()
    return
  end if
end if
if (Verbose) then
  call CWTime(CPUEnr2,WallEnr2)
  call Cho_PrtTim('Cholesky MP2 energy',CPUEnr2,CPUEnr1,WallEnr2,WallEnr1,iFmt)
end if

call Finish_this()

return

contains

! Exit.
! -----
subroutine Finish_this()

  real(kind=wp) :: Diff

  Diff = abs(Check(1)-Chk_Mem_ChoMP2)
  if (Diff > Tol) then
    write(u6,*) SecNam,': Memory Boundary Error!'
    if (irc == 0) irc = -9999
  end if
  if (Verbose) then
    call CWTime(CPUTot2,WallTot2)
    call Cho_PrtTim('Cholesky MP2',CPUTot2,CPUTot1,WallTot2,WallTot1,iFmt)
  end if
  call ChoMP2_deallocate(irc)
  call mma_deallocate(Check)

end subroutine Finish_this

end subroutine ChoMP2_Drv
