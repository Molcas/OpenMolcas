************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Aces_Gamma()
      use Basis_Info, only: nBas
      use Aces_Stuff
      use Index_arrays, only: iSO2Sh
      use Real_Info, only: CutInt
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "setup.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
       Integer, Allocatable:: iTable(:,:)
       Real*8, Allocatable:: Buf(:), Bin3(:,:,:)
*                                                                      *
************************************************************************
*                                                                      *
      nShell=mSkal
      nPair=nShell*(nShell+1)/2
      nQuad=nPair*(nPair+1)/2
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate Table Of Content for half sorted gammas.
*
      Call mma_allocate(G_Toc,nQuad,Label='G_Toc')
*
*---- Tabel SO to contigeous index
*
      Call mma_allocate(SO2cI,2,nSOs,Label='SO2cI')
*
*     Both should be deallocated in CloseP!
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate table with information regarding the symmetry blocks
*     of the gammas as stored in Aces 2 format
*
      If (nIrrep.eq.8) nBlocks=106
      If (nIrrep.eq.4) nBlocks= 19
      If (nIrrep.eq.2) nBlocks=  4
      If (nIrrep.eq.1) nBlocks=  1
      Call mma_Allocate(iTable,6,nBlocks,Label='iTable')
      Call Gamma_Blocks(iTable,nBlocks,nIrrep)
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate memory for read buffer
*
      Call GetMem('Buff','Max','Real',iDummy,MaxMem)
      nReq=0
      Do iBlock = 1, nBlocks
         iType=iTable(1,iBlock)
         iIrrep_A=iTable(2,iBlock)
         iIrrep_B=iTable(3,iBlock)
         iIrrep_C=iTable(4,iBlock)
         iIrrep_D=iTable(5,iBlock)
*
         nA=nBas(iIrrep_A)
         nB=nBas(iIrrep_B)
         nC=nBas(iIrrep_C)
         nD=nBas(iIrrep_D)
         If (iType.eq.1 .or. iType.eq.2) Then
            nAB = nA*(nA+1)/2
            nCD = nC*(nC+1)/2
         Else
            nAB=nA*nB
            nCD=nC*nD
         End If
*
         nReq=Max(nReq,nAB*nCD)
      End Do
      nReq=Min(MaxMem/4,nReq)
      Call mma_allocate(Buf,nReq,Label='Buf')
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate bins for shell quadruplets
*
      Call GetMem('Buff','Max','Real',iDummy,MaxMem)
      lBin=Min(MaxMem/(2*nQuad),1024)
      Call mma_allocate(Bin3,2,lBin,nQuad,Label='Bin3')
*                                                                      *
************************************************************************
*                                                                      *
*---- Open the bin file with half sorted gammas.
*
      LuGamma=60
      LuGamma=isfreeunit(LuGamma)
      Call DaName_MF(LuGamma,'GAMMA')
*                                                                      *
************************************************************************
*                                                                      *
*---- Read the blocks off the Aces 2 file and put into half sorted bin
*     file. The second half sort is done on the fly as needed.
*
      Call Read_Blocks(iTable,nBlocks,nBas,nIrrep,Buf,nReq,
     &                 iSO2Sh,nSOs,Bin3,lBin,nQuad,G_Toc,SO2cI,CutInt)
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory
*
      Call mma_deallocate(Bin3)
      Call mma_deallocate(Buf)
      Call mma_deallocate(iTable)
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate buffer for reading the bins
*
      Call mma_allocate(Bin,2,lBin,Label='Bin')
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
