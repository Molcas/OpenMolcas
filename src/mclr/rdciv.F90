!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine RdCIV()
!***********************************************************************
!                                                                      *
!     Read the contents of the JOBIPH file.                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use gugx, only: SGS, CIS, EXS
use MCLR_Data, only: ERAS, LuCIV
use MCLR_Data, only: LuJob
use input_mclr, only: nSym, lRoots, nCSF, nConf, nRS1, nRS2, nRS3, State_Sym, iSpin, iTOC, nActEl, nElec3, nHole1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero

implicit none
#include "rasdim.fh"
real*8, allocatable :: OCIvec(:), Tmp(:)
integer iDisk, iDisk1, i, Iter
real*8 Temp

call DaName(LuCIV,'ROOTS')
iDisk = 0
!----------------------------------------------------------------------*
!     Load the CI vector for all roots                                 *
!----------------------------------------------------------------------*

call mma_Allocate(OCIvec,nConf,Label='OCIvec')

iDisk = iToc(4)
idisk1 = 0
do i=1,lroots
  call dDaFile(LuJob,2,OCIvec,nConf,iDisk)
  call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,OCIvec,1,State_Sym,State_Sym)
  NCSF(1:nSym) = CIS%NCSF(1:nSym)
  NCONF = CIS%NCSF(State_Sym)
  call mkGuga_Free(SGS,CIS,EXS)

  call dDafile(LuCIV,1,OCIvec,nconf,iDisk1)
end do

call mma_deAllocate(OCIvec)
!----------------------------------------------------------------------*
!     Load state energy                                                *
!----------------------------------------------------------------------*

call mma_allocate(Tmp,mxRoot*mxIter,Label='Tmp')
iDisk = iToc(6)
call dDaFile(LuJob,2,Tmp,mxRoot*mxIter,iDisk)
do i=0,lroots-1
  ERAS(i+1) = Zero
  do iter=1,mxIter
    Temp = Tmp(iter*mxRoot+i)
    if (Temp /= Zero) ERAS(i+1) = Temp
  end do
end do
call mma_deallocate(Tmp)
!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*

end subroutine RdCIV
