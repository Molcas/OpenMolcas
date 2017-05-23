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
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "shinf.fh"
#include "setup.fh"
#include "WrkSpc.fh"
#include "aces_gamma.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iTable(i,j)=iWork(ipiTab-1+(j-1)*6+i)
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
      Call GetMem('G_Toc','Allo','Real',ipG_Toc,nQuad)
*
*---- Tabel SO to contigeous index
*
      Call GetMem('SO2cI','Allo','Inte',ipSO2cI,2*nSOs)
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
      Call GetMem('iTable','Allo','Inte',ipiTab,6*nBlocks)
      Call Gamma_Blocks(iWork(ipiTab),nBlocks,nIrrep)
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
      Call GetMem('Buff','Allo','Real',ipBuf,nReq)
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate bins for shell quadruplets
*
      Call GetMem('Buff','Max','Real',iDummy,MaxMem)
      lBin=Min(MaxMem/(2*nQuad),1024)
*     Write (*,*) 'lBin=',lBin
      Call GetMem('Bin','Allo','Real',ipBin,2*lBin*nQuad)
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
      Call Read_Blocks(iWork(ipiTab),nBlocks,nBas,nIrrep,
     &                 iOffSO,Work(ipBuf),nReq,
     &                 iWork(ipSOSh),nSOs,Work(ipBin),lBin,nQuad,
     &                 Work(ipG_Toc),iWork(ipSO2cI),CutInt)
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory
*
      Call GetMem('Bin','Free','Real',ipBin,2*lBin*nQuad)
      Call GetMem('Buff','Free','Real',ipBuf,nReq)
      Call GetMem('iTable','Free','Inte',ipiTab,6*nBlocks)
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate buffer for reading the bins
*
      Call GetMem('Bin','Allo','Real',ipBin,2*lBin)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
