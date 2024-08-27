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
use iSD_data, only: iSD, nShBf, iShOff, nShBfMx, iCntr, iSh2Sh, nShIrp, iSO2Sh
use Basis_Info, only: nBas, nBas_Aux
use stdalloc, only: mma_allocate
use BasisMode, only: Basis_Mode, Auxiliary_Mode

implicit none
integer nSkal, nSym, nSOs
integer nShOff(0:7), iTmp(1), iSkal, iAO, iCmp, iSO, iPtr, i, iRP, nShBfi, iSOB

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
      write(6,*) nShBF(irp,iskal)
      write(6,*) nShOff(irp)
      write(6,*) iShOff(irp,iSkal)
      write(6,*) 'PROGRAMMING ERROR IN SHELL_SIZES: SHELLS NOT CONTIGUOUS. IRP=',irp,'  ISKAL=',iskal
      call Abend()
    end if
#   endif
    iShOff(irp,iSkal) = nShOff(irp)
    nShOff(irp) = nShOff(irp)+nShBF(irp,iSkal)
  end do
# ifdef _DEBUGPRINT_
  write(6,'(A)') 'nShBF'
  write(6,'(8I4)') iSkal,(nShBF(irp,iSkal),irp=0,nSym-1)
  write(6,'(A)') 'iShOff'
  write(6,'(8I4)') iSkal,(iShOff(irp,iSkal),irp=0,nSym-1)
# endif
end do

! and now set up SO-Shell and Shell-Psudoshell index vectors...

iTmp = 0 ! Use iTmp to get round compiler bug on some machines
call ICopy(nSym,iTmp,0,nShIrp,1)

iTmp = -9999999
call ICopy(nSOs,iTmp,0,iSO2Sh,1)
call ICopy(nSkal*nSym,iTmp,0,iSh2Sh,1)

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
write(6,'(A,I4)') 'max shell size:',nShBFMx
write(6,'(A)') '# of shells contributing to each irrep:'
write(6,'(8I4)') (nShIrp(irp),irp=0,nSym-1)
write(6,'(A)') '# shell-psudoshell map vector:'
do irp=0,nSym-1
  write(6,'(A4,2X,I4,2X,A4,2X,8I4)') 'irp=',irp,'map:',(iSh2Sh(irp,iSkal),iSkal=1,nSkal)
end do
write(6,'(A)') 'SO-index to shell map vector:'
write(6,'(50I4)') (iSO2Sh(iSO),iSO=1,nSOs)
#endif
return

end subroutine SOFSh1
