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
! Copyright (C) 2007,2008, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine ChoMP2_VectorMO2AO(iTyp,Delete,BaseName_AO,CMO,DoDiag,Diag,lDiag,lU_AO,irc)
!
! Thomas Bondo Pedersen, Dec. 2007 - Jan. 2008.
!
! Purpose: backtransform vectors from MO to AO basis.
!
! Input:
!    iTyp......... specifies type of vectors (iTyp is used to open
!                  MO vector files through ChoMP2_OpenF()).
!    Delete....... Flag specifiyng whether MO vector files are to
!                  be deleted before exiting this routine.
!    BaseName_AO.. Base name (3 characters) for the files containing
!                  AO vectors. Symmetry index will be appended!
!    CMO.......... MO coefficient array
!    DoDiag....... if .True., calculate AO diagonal elements as
!                  D(ab) = sum_J L(J,ab)*L(J,ab)
! Output:
!    Diag......... Contains the diagonal if requested (flag DoDiag).
!    lDiag........ Dimension of Diag
!    lU_AO........ Array containing units of open AO vector files.
!    irc.......... return code.
!                  = 0 if successful, non-zero otherwise
!                  (must be checked by caller).
!
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!*** CHOLESKY INFORMATION MUST BE INITIALIZED WHEN CALLING THIS ROUTINE.
!    --> I.e. Cho_X_Init() must have been called.
!    --> It is actually sufficient that nBas + nSym in
!        cholesky are available (which they are when Cho_X_Init()
!        has been called).
!*** CHOLESKY MP2 INFORMATION MUST BE INITIALIZED WHEN CALLING THIS
!    ROUTINE.
!    --> I.e. ChoMP2_Setup() must have been called.
!*** AO VECTORS ARE STORED IN LOWER TRIANGULAR [M(I,J), I >= J] FORMAT.
!    --> I.e. vectors are stored as L(J,ab) where a>=b
!*** DIAGONAL IS STORED IN LOWER TRIANGULAR FORMAT.
!    --> I.e. diagonal is stored as D(ab) where a>=b (same as vectors).
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************

use Symmetry_Info, only: Mul
use Cholesky, only: nBas, nSym
use ChoMP2, only: nAOVir, nT1AOT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iTyp, lDiag
logical(kind=iwp), intent(in) :: Delete, DoDiag
character(len=3), intent(in) :: BaseName_AO
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(out) :: Diag(lDiag)
integer(kind=iwp), intent(out) :: lU_AO(*), irc
integer(kind=iwp) :: iClose, iCount, iOpen, iSym, iSyma, iSymb
character(len=4) :: FullName_AO
real(kind=wp), allocatable :: COcc(:), CVir(:)
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: Debug = _DBG_
character(len=*), parameter :: SecNam = 'ChoMP2_VectorMO2AO'

! Initializations.
! ----------------

irc = 0
lU_AO(1:nSym) = -999999
if (DoDiag) then
  iCount = 0
  do iSym=1,nSym
    do iSymb=1,nSym
      iSyma = Mul(iSymb,iSym)
      iCount = iCount+nBas(iSyma)*nBas(iSymb)
    end do
  end do
  if (iCount /= lDiag) then
    write(u6,*) SecNam,': WARNING: inconsistent diagonal allocation!'
    if (iCount > lDiag) then
      write(u6,*) '   - insufficient memory, will return now...'
      irc = 1
      return
    else
      write(u6,*) '   - sufficient memory, going to continue...'
    end if
  end if
end if

! Reorder CMO. This also removes frozen orbitals.
! -----------------------------------------------

call mma_allocate(COcc,nT1AOT(1),Label='COcc')
call mma_allocate(CVir,nAOVir(1),Label='CVir')
call ChoMP2_MOReOrd(CMO,COcc,CVir)

! Backtransform.
! --------------

call ChoMP2_BackTra(iTyp,COcc,CVir,BaseName_AO,DoDiag,Diag)

! Open AO vector files (i.e. get units to return).
! ------------------------------------------------

do iSym=1,nSym
  write(FullName_AO,'(A3,I1)') BaseName_AO,iSym
  lU_AO(iSym) = 7
  call daName_MF_WA(lU_AO(iSym),FullName_AO)
end do

! Debug: check backtransformation.
! --------------------------------

if (Debug) call ChoMP2_CheckBackTra(iTyp,COcc,CVir,lU_AO)

! Delete MO files if requested.
! -----------------------------

if (Delete) then
  iOpen = 1
  iClose = 3
  do iSym=1,nSym
    call ChoMP2_OpenF(iOpen,iTyp,iSym)
    call ChoMP2_OpenF(iClose,iTyp,iSym)
  end do
end if

! Deallocate and exit.
! --------------------

call mma_deallocate(CVir)
call mma_deallocate(COcc)

end subroutine ChoMP2_VectorMO2AO
