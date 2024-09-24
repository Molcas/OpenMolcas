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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine MkEorb_Inner(FockAO,nFck,CMO,nCMO,Eorb,nEorb,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
!  This routine calculates the diagonal elements of the MO Fock matrix *
!  (orbital energies).                                                 *
!                                                                      *
!  Input:                                                              *
!    FockAO  Fock matrix in AO basis                                   *
!    CMO     Orbitals                                                  *
!                                                                      *
!  Output:                                                             *
!    Eorb    Orbital energies.                                         *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero

implicit none
!----------------------------------------------------------------------*
! Dummy arguments.                                                     *
!----------------------------------------------------------------------*
integer nFck, nCMO, nEOrb
real*8 FockAO(nFck)
real*8 CMO(nCMO)
real*8 EOrb(nEOrb)
integer nSym
integer nBas(nSym)
integer nOrb(nSym)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
real*8 t
real*8, dimension(:), allocatable :: FckSqr
integer iSym
integer iBas
integer jBas
integer iOrb
integer MaxTri
integer MaxSqr
integer iOffTri
integer iOffCMO
integer npFckSqr
integer indE
integer indF
integer indx
integer jndx

!----------------------------------------------------------------------*
! Some preliminary setup.                                              *
!----------------------------------------------------------------------*
MaxTri = 0
MaxSqr = 0
do iSym=1,nSym
  MaxTri = max(MaxTri,nBas(iSym)*(nBas(iSym)+1)/2)
  MaxSqr = max(MaxSqr,nBas(iSym)*nBas(iSym))
end do
!----------------------------------------------------------------------*
! Allocate matrices.                                                   *
!----------------------------------------------------------------------*
npFckSqr = MaxSqr
call mma_allocate(FckSqr,npFckSqr,Label='FckSqr')
!----------------------------------------------------------------------*
! Do compute orbital energies                                          *
!----------------------------------------------------------------------*
iOffTri = 0
iOffCMO = 0
indE = 1
do iSym=1,nSym
  if (nOrb(iSym) > 0) then
    call Square(FockAO(1+iOffTri),FckSqr,1,nBas(iSym),nBas(iSym))
    do iOrb=1,nOrb(iSym)
      t = Zero
      indF = 1
      do iBas=1,nBas(iSym)
        do jBas=1,nBas(iSym)
          indx = iBas+(iOrb-1)*nBas(iSym)+iOffCMO
          jndx = jBas+(iOrb-1)*nBas(iSym)+iOffCMO
          t = t+CMO(indx)*CMO(jndx)*FckSqr(indF)
          indF = indF+1
        end do
      end do
      EOrb(indE) = t
      indE = indE+1
    end do
  end if
  iOffTri = iOffTri+nBas(iSym)*(nBas(iSym)+1)/2
  iOffCMO = iOffCMO+nBas(iSym)*nOrb(iSym)
end do
!----------------------------------------------------------------------*
! Debug print                                                          *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
indE = 1
do iSym=1,nSym
  if (nOrb(iSym) > 0) call RecPrt('mkeor: Orbital energies','(20F10.4)',Eorb(indE),1,nOrb(iSym))
  indE = indE+nOrb(iSym)
end do
#endif
!----------------------------------------------------------------------*
! Deallocate matrices.                                                 *
!----------------------------------------------------------------------*
call mma_deallocate(FckSqr)
!----------------------------------------------------------------------*
! Done.                                                                *
!----------------------------------------------------------------------*
return

end subroutine MkEorb_Inner
