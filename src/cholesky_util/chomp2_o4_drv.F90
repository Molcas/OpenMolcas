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
! Copyright (C) 2008, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_O4_Drv(irc,EMP2,CMO,EOcc,EVir)
!
! Thomas Bondo Pedersen, Jan. 2008.
! - based on ChoMP2_Drv() by T. B. Pedersen.
!
! Purpose: driver for computing the MP2 energy correction EMP2
!          using Cholesky decomposed two-electron integrals
!          in a quartic scaling fashion.
!          Input must have been processed and MO coefficients
!          and orbital energies must be passed as arguments.
!
! Notes:
!
!   - all MO Cholesky vector files generated here are deleted before
!     exit, except for error terminations (i.e. no cleanup actions
!     are taken!)

use Symmetry_Info, only: Mul
use Cholesky, only: nBas, nSym
use ChoMP2, only: iOcc, iT1am, iVir, nOcc, nT1am, nVir, Verbose
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(in) :: CMO(*), EOcc(*), EVir(*)
integer(kind=iwp) :: a, ai, i, iSym, iSyma, iSymb, iSymi, iTyp, kD0, kD1, kD2, lDiag, lU_AO(8)
real(kind=wp) :: CPUBT1, CPUBT2, CPUDec1, CPUDec2, CPUIni1, CPUIni2, CPUTot1, CPUTot2, CPUTra1, CPUTra2, DE, Ei, FracMem, WallBT1, &
                 WallBT2, WallDec1, WallDec2, WallIni1, WallIni2, WallTot1, WallTot2, WallTra1, WallTra2
logical(kind=iwp) :: Delete, DoAmpDiag
character(len=3) :: BaseName_AO
real(kind=wp), allocatable :: Check(:), Diag(:)
integer(kind=iwp), parameter :: iFmt = 0
real(kind=wp), parameter :: Chk_Mem_ChoMP2 = 0.123456789_wp, Tol = 1.0e-15_wp
logical(kind=iwp), parameter :: Delete_def = .true.
character(len=*), parameter :: SecNam = 'ChoMP2_O4_Drv'

#ifdef _DEBUGPRINT_
Verbose = .true.
#endif
if (Verbose) call CWTime(CPUTot1,WallTot1)

! Initializations.
! ----------------

irc = 0

EMP2 = Zero

if (Verbose) call CWTime(CPUIni1,WallIni1)

call mma_allocate(Check,1,Label='Check')
Check(1) = Chk_Mem_ChoMP2

FracMem = Zero ! no buffer allocated
call Cho_X_Init(irc,FracMem)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_X_Init returned ',irc
  call SysAbendMsg(SecNam,'Cholesky initialization error',' ')
end if

!-TBP:
! Frankie,
! The setup is still the same here (i.e. batching etc. is included)
! I don't use the batching info at all, though, so it should be save
! to remove it - unless you need it, of course.

call ChoMP2_Setup(irc)
if (irc /= 0) then
  write(u6,*) SecNam,': ChoMP2_Setup returned ',irc
  call Finish_this()
  return
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
! Compute also amplitude diagonal here.
! ----------------------------------------------------------

if (Verbose) call CWTime(CPUTra1,WallTra1)

lDiag = sum(nT1am(1:nSym))

call mma_allocate(Diag,lDiag,Label='Diag')

call ChoMP2_TraDrv(irc,CMO,Diag,.true.)
if (irc /= 0) then
  write(u6,*) SecNam,': ChoMP2_TraDrv returned ',irc
  call Finish_this()
  return
end if
kD0 = 0
do iSym=1,nSym
  do iSymi=1,nSym
    iSyma = Mul(iSymi,iSym)
    kD1 = kD0+iT1Am(iSyma,iSymi)
    do i=1,nOcc(iSymi)
      kD2 = kD1+nVir(iSyma)*(i-1)
      Ei = EOcc(iOcc(iSymi)+i)
      do a=1,nVir(iSyma)
        ai = kD2+a
        DE = Two*(EVir(iVir(iSyma)+a)-Ei)
        Diag(ai) = Diag(ai)/DE
      end do
    end do
  end do
  kD0 = kD0+nT1Am(iSym)
end do

if (Verbose) then
  call CWTime(CPUTra2,WallTra2)
  call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,CPUTra1,WallTra2,WallTra1,iFmt)
end if

! Decompose MP2 amplitudes (times -1).
! ------------------------------------

if (Verbose) call CWTime(CPUDec1,WallDec1)

!-TBP:
! Frankie,
! I just modified the decomposition slightly so that it treats
! amplitudes, too. You should be aware, though, that result vectors
! are always written on the same files, so if you do another
! decomposition of, say, integrals (or squared integrals), then
! the vector files will be overwritten!!
! The number of vectors is always written to nMP2Vec(iSym) in
! ChoMP2 - this is overwritten too, if you do another CD!!

Delete = Delete_def ! delete transf. vector files after dec.
call ChoMP2_DecDrv(irc,Delete,Diag,'Amplitudes')
if (irc /= 0) then
  write(u6,*) SecNam,': ChoMP2_DecDrv returned ',irc
  call SysAbendMsg(SecNam,'MP2 decomposition failed!',' ')
end if
call mma_deallocate(Diag)

if (Verbose) then
  call CWTime(CPUDec2,WallDec2)
  call Cho_PrtTim('Cholesky MP2 decomposition',CPUDec2,CPUDec1,WallDec2,WallDec1,iFmt)
end if

! Backtransform amplitude vectors to AO basis.
! Calculate also backtransformed amplitude diagonal.
! --------------------------------------------------

if (Verbose) call CWTime(CPUBT1,WallBT1)

iTyp = 2 ! type of MO vectors (i.e. amp. vectors)
Delete = Delete_def ! delete amp. vectors after backtransf.
BaseName_AO = 'AAO' ! basename for multifiles of AO vectors
DoAmpDiag = .true. ! calculate backtransf. amp. diagonal
lDiag = 0
do iSym=1,nSym
  do iSymb=1,nSym
    iSyma = Mul(iSymb,iSym)
    lDiag = lDiag+nBas(iSyma)*nBas(iSymb)
  end do
end do

call mma_allocate(Diag,lDiag,Label='Diag')
call ChoMP2_VectorMO2AO(iTyp,Delete,BaseName_AO,CMO,DoAmpDiag,Diag,lDiag,lU_AO,irc)
if (irc /= 0) then
  write(u6,*) SecNam,': ChoMP2_VectorMO2AO returned ',irc
  call SysAbendMsg(SecNam,'MP2 amplitude vector backtransformation failed!',' ')
end if
call mma_deallocate(Diag)

if (Verbose) then
  call CWTime(CPUBT2,WallBT2)
  call Cho_PrtTim('Cholesky MP2 backtransformation',CPUBT2,CPUBT1,WallBT2,WallBT1,iFmt)
end if

!-TBP:
! Frankie,
! You can now read the AO vectors from the units lU_AO(iSym)
! using ddaFile(). The files are word-addressable so that it
! is possible to read from the file as if addressing an array.
! The number of vectors is nMP2Vec(iSym) stored in ChoMP2.
! To save memory, you may want to finalize Cholesky info before
! continuing with the backtransformed vectors - but remember
! that all Cholesky information (from the AO integral CD) is
! then lost (f.ex. NumCho(iSym) becomes useless).

! Finalize Cholesky info.
! -----------------------

call Cho_X_Final(irc)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_X_Final returned ',irc
  call Finish_this()
  return
end if

! Close and delete files containing backtransformed amplitude vectors.
! --------------------------------------------------------------------

do iSym=1,nSym
  call daEras(lU_AO(iSym))
end do

call Finish_this()

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

  call mma_deallocate(Check)

end subroutine Finish_this

end subroutine ChoMP2_O4_Drv
