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
! Copyright (C) 1995, Martin Schuetz                                   *
!***********************************************************************

!#define _CHECK_
subroutine SOFSh1(nSkal,nSym,nSOs)
!***********************************************************************
! This Module contains subroutines which are used to compute info on   *
! the size of the SO integral symmetry blocks for direct integral      *
! transformation                                                       *
!                                                                      *
! SubRoutine SOFSh1                                                    *
!  -> compute (1) # SO functions in irrep for all shells iShell        *
!             (2) position of 1st component of shell in irrep for      *
!                 all shells in all irreps                             *
!             (3) map vector between shells and psedoshells for each   *
!                 irrep (not any shell contributes to all irreps)      *
!             (4) map vector between SO indices and shells             *
! relevant data is declared and passed in common "inftra.fh"           *
!----------------------------------------------------------------------*
!     written by:                                                      *
!     M. Schuetz                                                       *
!     University of Lund, Sweden, 1995                                 *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use iSD_data, only: iCntr, iSD, iSh2Sh, iShOff, iSO2Sh, nShBf, nShBfMx, nShIrp
use Basis_Info, only: nBas, nBas_Aux
use BasisMode, only: Auxiliary_Mode, Basis_Mode
use stdalloc, only: mma_allocate
use Definitions, only: iwp
#if defined (_DEBUGPRINT_) || defined (_CHECK_)
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nSkal, nSym, nSOs
integer(kind=iwp) :: i, iAO, iCmp, iPtr, iRP, iSkal, iSO, iSOB, nShBfi, nShOff(0:7)

! Allocate all memory

call mma_allocate(nShBF,[0,nSym-1],[1,nSkal],Label='nShBF')
call mma_allocate(iShOff,[0,nSym-1],[1,nSkal],Label='iShOff')
call mma_allocate(iSh2Sh,[0,nSym-1],[1,nSkal],Label='iSh2Sh')
call mma_allocate(iSO2Sh,nSOs,Label='iSO2Sh')
call mma_allocate(iCntr,nSkal,Label='iCntr')

! Initialize
nShBF(:,:) = 0
iShOff(:,:) = 9999999
nShOff(:) = 1

do iSkal=1,nSkal
  iAO = iSD(7,iSkal)
  iCmp = iSD(2,iSkal)
  icntr(iSkal) = iSD(10,iSkal)

  ! loop over components of shell...

  do i=1,iCmp
    ! loop over irreps...
    do irp=0,nSym-1
      if (iAOtSO(iAO+i,irp) > 0) then
        nShBF(irp,iSkal) = nShBF(irp,iSkal)+iSD(3,iSkal)
#       ifdef _CHECK_
        if (Basis_Mode == Auxiliary_Mode) then
          iShOff(irp,iSkal) = min(iShOff(irp,iSkal),iAOtSO(iAO+i,irp)-nBas(irp))
        else
          iShOff(irp,iSkal) = min(iShOff(irp,iSkal),iAOtSO(iAO+i,irp))
        end if
#       endif
      end if
    end do
  end do
  do irp=0,nSym-1
#   ifdef _CHECK_
    if ((nShBF(irp,iskal) /= 0) .and. (nShOff(irp) /= iShOff(irp,iSkal))) then
      call WarningMessage(2,'PROGRAMMING ERROR IN SHELL_SIZES')
      write(u6,*) nShBF(irp,iskal)
      write(u6,*) nShOff(irp)
      write(u6,*) iShOff(irp,iSkal)
      write(u6,*) 'PROGRAMMING ERROR IN SHELL_SIZES: SHELLS NOT CONTIGUOUS. IRP=',irp,'  ISKAL=',iskal
      call Abend()
    end if
#   endif
    iShOff(irp,iSkal) = nShOff(irp)
    nShOff(irp) = nShOff(irp)+nShBF(irp,iSkal)
  end do
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'nShBF'
  write(u6,'(8I4)') iSkal,(nShBF(irp,iSkal),irp=0,nSym-1)
  write(u6,'(A)') 'iShOff'
  write(u6,'(8I4)') iSkal,(iShOff(irp,iSkal),irp=0,nSym-1)
# endif
end do

! and now set up SO-Shell and Shell-Psudoshell index vectors...

nShIrp(0:nSym-1) = 0
iSO2Sh(:) = -9999999
iSh2Sh(:,:) = -9999999

! Loop over irreps...

iptr = 0
nShBFMx = 0
do irp=0,nSym-1
  do iSkal=1,nSkal

    nShBFi = nShBF(irp,iSkal)
    nShBFMx = max(nShBFMx,nShBFi)
    iSOb = iShOff(irp,iSkal)
    do iSO=iSOb,iSOb+nShBFi-1
      if (iSO > nSOs) then
        call WarningMessage(2,' Fucked again!')
        call Quit_OnUserError()
      end if
      iSO2Sh(iptr+iSO) = iSkal
    end do

    if (nShBFi > 0) then
      nShIrp(irp) = nShIrp(irp)+1
      iSh2Sh(irp,iSkal) = nShIrp(irp)
    end if
  end do
  if (Basis_Mode == Auxiliary_Mode) then
    iptr = iptr+nBas_Aux(irp)
  else
    iptr = iptr+nBas(irp)
  end if
end do
#ifdef _DEBUGPRINT_
write(u6,'(A,I4)') 'max shell size:',nShBFMx
write(u6,'(A)') '# of shells contributing to each irrep:'
write(u6,'(8I4)') (nShIrp(irp),irp=0,nSym-1)
write(u6,'(A)') '# shell-psudoshell map vector:'
do irp=0,nSym-1
  write(u6,'(A4,2X,I4,2X,A4,2X,8I4)') 'irp=',irp,'map:',(iSh2Sh(irp,iSkal),iSkal=1,nSkal)
end do
write(u6,'(A)') 'SO-index to shell map vector:'
write(u6,'(50I4)') (iSO2Sh(iSO),iSO=1,nSOs)
#endif
return

end subroutine SOFSh1
