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

subroutine Aces_Gamma()

use Index_Functions, only: nTri_Elem
use setup, only: mSkal, nSOS
use Basis_Info, only: nBas
use PSO_Stuff, only: Bin, G_ToC, lBin, LuGamma, SO2CI
use iSD_data, only: iSO2Sh
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iBlock, iIrrep_A, iIrrep_B, iIrrep_C, iIrrep_D, iType, MaxMem, nA, nAB, nB, nBlocks, nC, nCD, nD, nPair, &
                     nQUad, nReq, nShell
integer(kind=iwp), allocatable :: iTable(:,:)
real(kind=wp), allocatable :: Buf(:), Bin3(:,:,:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
nShell = mSkal
nPair = nTri_Elem(nShell)
nQuad = nTri_Elem(nPair)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate Table Of Content for half sorted gammas.

call mma_allocate(G_Toc,nQuad,Label='G_Toc')

! Table SO to contiguous index

call mma_allocate(SO2cI,2,nSOs,Label='SO2cI')

! Both should be deallocated in CloseP!
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate table with information regarding the symmetry blocks
! of the gammas as stored in Aces 2 format

if (nIrrep == 8) nBlocks = 106
if (nIrrep == 4) nBlocks = 19
if (nIrrep == 2) nBlocks = 4
if (nIrrep == 1) nBlocks = 1
call mma_Allocate(iTable,6,nBlocks,Label='iTable')
call Gamma_Blocks(iTable,nBlocks,nIrrep)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for read buffer

call mma_MaxDBLE(MaxMem)
nReq = 0
do iBlock=1,nBlocks
  iType = iTable(1,iBlock)
  iIrrep_A = iTable(2,iBlock)
  iIrrep_B = iTable(3,iBlock)
  iIrrep_C = iTable(4,iBlock)
  iIrrep_D = iTable(5,iBlock)

  nA = nBas(iIrrep_A)
  nB = nBas(iIrrep_B)
  nC = nBas(iIrrep_C)
  nD = nBas(iIrrep_D)
  if ((iType == 1) .or. (iType == 2)) then
    nAB = nTri_Elem(nA)
    nCD = nTri_Elem(nC)
  else
    nAB = nA*nB
    nCD = nC*nD
  end if

  nReq = max(nReq,nAB*nCD)
end do
nReq = min(MaxMem/4,nReq)
call mma_allocate(Buf,nReq,Label='Buf')
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate bins for shell quadruplets

call mma_MaxDBLE(MaxMem)
lBin = min(MaxMem/(2*nQuad),1024)
call mma_allocate(Bin3,2,lBin,nQuad,Label='Bin3')
!                                                                      *
!***********************************************************************
!                                                                      *
! Open the bin file with half sorted gammas.

LuGamma = 60
LuGamma = isfreeunit(LuGamma)
call DaName_MF(LuGamma,'GAMMA')
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the blocks off the Aces 2 file and put into half sorted bin
! file. The second half sort is done on the fly as needed.

call Read_Blocks(iTable,nBlocks,nBas,nIrrep,Buf,nReq,iSO2Sh,nSOs,Bin3,lBin,nQuad,G_Toc,SO2cI,CutInt)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory

call mma_deallocate(Bin3)
call mma_deallocate(Buf)
call mma_deallocate(iTable)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate buffer for reading the bins

call mma_allocate(Bin,2,lBin,Label='Bin')
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Aces_Gamma
