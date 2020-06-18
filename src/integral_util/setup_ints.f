************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Roland Lindh                                     *
************************************************************************
      Subroutine SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
************************************************************************
*                                                                      *
*     Object: to set up data and allocate memory for integral calcula- *
*             tions. The whole data structure is hidden to the user.   *
*                                                                      *
*     nSkal(output): returns the number of shells                      *
*     Indexation(input): logical flag to initiate index tables         *
*     ThrAO(input): if ThrAO.ne.Zero CutInt is reset to ThrAO          *
*                                                                      *
*     Author: Roland Lindh, Chemical Physics, University of Lund,      *
*             Sweden. January '98.                                     *
************************************************************************
      use Her_RW
      use vRys_RW
      use iSD_data
      use k2_arrays
      use index_arrays
      Implicit Real*8 (a-h,o-z)
      External CmpctR, CmpctS
#include "itmax.fh"
#include "info.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "lundio.fh"
#include "setup.fh"
#include "real.fh"
#include "shinf.fh"
#include "status.fh"
#include "ndarray.fh"
*
      Logical DoFock, DoGrad, Indexation
*
      If (ERI_Status.eq.Active) Then
        Call Nr_Shells(nSkal)
        Return
      End If
      ERI_Status=Active
*     Call QEnter('S_I')
*                                                                      *
************************************************************************
*                                                                      *
      if(thrao.ne.Zero) CutInt=ThrAO
*
*.....Compute the total number of symmetry adapted basis functions
*
      nSOs = 0
      Do iIrrep = 0, nIrrep-1
         If (Basis_Mode.eq.Valence_Mode) Then
            nSOs = nSOs + nBas(iIrrep)
         Else If (Basis_Mode.eq.Auxiliary_Mode) Then
            nSOs = nSOs + nBas_Aux(iIrrep)
         Else If (Basis_Mode.eq.With_Auxiliary_Mode) Then
            nSOs = nSOs + nBas(iIrrep) + nBas_Aux(iIrrep)
         End If
      End Do
*
*.....Generate a two-dimensional array of the length of nSOs.
*     The first entry gives the irrep of a SO and the second entry
*     gives the relative index of a SO in its irrep.
*
      Call mma_allocate(iSOSym,2,nSOs,Label='iSOSym')
      iSOs = 1
      nBas_iIrrep=0
      Do iIrrep = 0, nIrrep-1
         If (Basis_Mode.eq.Valence_Mode) Then
            nBas_iIrrep=nBas(iIrrep)
         Else If (Basis_Mode.eq.Auxiliary_Mode) Then
            nBas_iIrrep=nBas_Aux(iIrrep)
         Else If (Basis_Mode.eq.With_Auxiliary_Mode) Then
            nBas_iIrrep=nBas(iIrrep)+nBas_Aux(iIrrep)
         End If
         Do i = 1, nBas_iIrrep
            iSOSym(1,iSOs)=iIrrep          ! Irreducible reps.
            iSOSym(2,iSOs)=i               ! Relative index in irrep.
            iSOs = iSOs + 1
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*.....Compute the number of shells and set up the shell information
*     tables(iSD).
*
      Call Nr_Shells(nSkal)
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     allocate Integer memory for resulting SO info...
*     memory basepointers are declared in inftra common block
*
      If (Indexation) Then
         Indexation_Status=Active
         Call GetMem('nShBF','ALLO','Inte',ipShBF,nSkal*nIrrep)
         Call GetMem('ShLwC','ALLO','Inte',ipShLC,nSkal*nIrrep)
         Call GetMem('ShPSh','ALLO','Inte',ipShSh,nSkal*nIrrep)
         Call mma_allocate(iSO2Sh,nSOs,Label='iSO2Sh')
         Call mma_allocate(iCntr,nSkal,Label='iCntr')
         Call SOFSh1(iWork(ipShBF),iWork(ipShLC),iWork(ipShSh),
     &               iSO2Sh,iCntr,nSkal,nIrrep,nSOs,nShIrp,nShBFmx)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate auxiliary array for symmetry transformation
*
      nAux = nIrrep**3
      If (Petite) nAux = 1
      Call mma_allocate(Aux,nAux,Label='Aux')
*                                                                      *
************************************************************************
*                                                                      *
*     Preallocate memory for k2 entities
*
      nZeta = MxPrm * MxPrm
      nEta  = MxPrm * MxPrm
      MemR=(nDArray-1)*nZeta + (nDArray-1)*nEta
      Call mma_allocate(Mem_DBLE,MemR,Label='Mem_DBLE')
      ipZeta=1
      MemI=nZeta+nEta+2
      Call mma_allocate(Mem_INT,MemI,Label='Mem_INT')
      ipiZet=1
*                                                                      *
************************************************************************
*                                                                      *
      If (DoFock) Then
         nFT=MxFT
      Else
         nFT=1
      End If
      Call mma_allocate(FT,MxFT,Label='FT')
*                                                                      *
************************************************************************
*                                                                      *
*     Precompute k2 entities
*
      If (lSchw) Then
         Call Drvk2(CmpctS,DoFock,DoGrad)
      Else
         Call Drvk2(CmpctR,DoFock,DoGrad)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call StatP(0)
      nUt=0
      iDisk=0
*                                                                      *
************************************************************************
*                                                                      *
*     Call QExit('S_I')
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Function iPD(iSO_,jSO_,iSOSym,nSOs)
#include "itmax.fh"
#include "info.fh"
      Integer iPD
      Integer iSOSym(2,nSOs)
*
      iPD = -999999
*
      iSO=Max(iSO_,jSO_)
      jSO=Min(iSO_,jSO_)
      iSym=iSOSym(1,iSO)
      iSOr=iSOSym(2,iSO)
      jSym=iSOSym(1,jSO)
      jSOr=iSOSym(2,jSO)
      If (iSym.eq.jSym) Then
          ij = iSOr*(iSOr-1)/2 + jSOr
      Else
          ij = (iSOr-1)*nBas(jSym) + jSOr
      End If
*
      iPD=ij
*
      Return
      End
