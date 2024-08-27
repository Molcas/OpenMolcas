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
! Copyright (C) 1998, Roland Lindh                                     *
!***********************************************************************

subroutine SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!***********************************************************************
!                                                                      *
!     Object: to set up data and allocate memory for integral calcula- *
!             tions. The whole data structure is hidden to the user.   *
!                                                                      *
!     nSkal(output): returns the number of shells                      *
!     Indexation(input): logical flag to initiate index tables         *
!     ThrAO(input): if ThrAO /= Zero CutInt is reset to ThrAO          *
!                                                                      *
!     Author: Roland Lindh, Chemical Physics, University of Lund,      *
!             Sweden. January '98.                                     *
!***********************************************************************

use setup, only: nSOs, nAux, MxPrm
use k2_arrays, only: nFT, MxFT, iSOSym, Aux, FT, create_braket_base
use Basis_Info, only: nBas, nBas_Aux
use Gateway_Info, only: CutInt, lSchw
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use stdalloc, only: mma_allocate
use BasisMode, only: Basis_Mode, Valence_Mode, Auxiliary_Mode, With_Auxiliary_Mode

implicit none
logical DoFock, DoGrad, Indexation
integer nSkal
real*8 ThrAO
external CmpctR, CmpctS
integer iIrrep, iSOs, nBas_iIrrep, i

if (allocated(iSOSym)) then
  call Nr_Shells(nSkal)
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (thrao /= Zero) CutInt = ThrAO

! Compute the total number of symmetry adapted basis functions

nSOs = 0
do iIrrep=0,nIrrep-1
  select case (Basis_Mode)
    case (Valence_Mode)
      nSOs = nSOs+nBas(iIrrep)
    case (Auxiliary_Mode)
      nSOs = nSOs+nBas_Aux(iIrrep)
    case (With_Auxiliary_Mode)
      nSOs = nSOs+nBas(iIrrep)+nBas_Aux(iIrrep)
  end select
end do

! Generate a two-dimensional array of the length of nSOs.
! The first entry gives the irrep of a SO and the second entry
! gives the relative index of a SO in its irrep.

call mma_allocate(iSOSym,2,nSOs,Label='iSOSym')
iSOs = 1
nBas_iIrrep = 0
do iIrrep=0,nIrrep-1
  select case (Basis_Mode)
    case (Valence_Mode)
      nBas_iIrrep = nBas(iIrrep)
    case (Auxiliary_Mode)
      nBas_iIrrep = nBas_Aux(iIrrep)
    case (With_Auxiliary_Mode)
      nBas_iIrrep = nBas(iIrrep)+nBas_Aux(iIrrep)
  end select
  do i=1,nBas_iIrrep
    iSOSym(1,iSOs) = iIrrep    ! Irreducible reps.
    iSOSym(2,iSOs) = i         ! Relative index in irrep.
    iSOs = iSOs+1
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the number of shells and set up the shell information
! tables(iSD).

call Nr_Shells(nSkal)
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
! allocate Integer memory for resulting SO info...
! memory basepointers are declared in inftra common block

if (Indexation) call SOFSh1(nSkal,nIrrep,nSOs)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate auxiliary array for symmetry transformation

nAux = nIrrep**3
if (nIrrep == 1) nAux = 1
call mma_allocate(Aux,nAux,Label='Aux')
!                                                                      *
!***********************************************************************
!                                                                      *
! Preallocate memory for k2 entities

call Create_BraKet_Base(MxPrm**2)
!                                                                      *
!***********************************************************************
!                                                                      *
if (DoFock) then
  nFT = MxFT
else
  nFT = 1
end if
call mma_allocate(FT,MxFT,Label='FT')
!                                                                      *
!***********************************************************************
!                                                                      *
! Precompute k2 entities

if (lSchw) then
  call Drvk2(CmpctS,DoFock,DoGrad)
else
  call Drvk2(CmpctR,DoFock,DoGrad)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine SetUp_Ints
