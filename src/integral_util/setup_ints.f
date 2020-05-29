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
      Implicit Real*8 (a-h,o-z)
      External CmpctR, CmpctS
#include "itmax.fh"
#include "info.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
#include "WrkSpc.fh"
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
      Call GetMem('iSOSym','Allo','Inte',ipiSOSym,nSOs*2)
      iSOs = ipiSOSym
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
            iWork(iSOs  )=iIrrep          ! Irreducible reps.
            iWork(iSOs+1)=i               ! Relative index in irrep.
            iSOs = iSOs + 2
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
         Call GetMem('SOShl','ALLO','Inte',ipSOSh,nSOs)
         Call GetMem('ICNTR','ALLO','Inte',ipicntr,nSkal)
         Call SOFSh1(iWork(ipShBF),iWork(ipShLC),iWork(ipShSh),
     &               iWork(ipSOSh),iWork(ipicntr),nSkal,nIrrep,nSOs,
     &               nShIrp,nShBFmx)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.Allocated(HerR) .or.
     &    .Not.Allocated(iHerR2)) Then
         Ind0_Status=Active
         Return
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate auxiliary array for symmetry transformation
*
      nAux = nIrrep**3
      If (Petite) nAux = 1
      Call GetMem('AuxBuf','ALLO','REAL',ipAux,nAux)
*                                                                      *
************************************************************************
*                                                                      *
*     Preallocate memory for k2 entities
*
      nZeta = MxPrm * MxPrm
      nEta  = MxPrm * MxPrm
      MemR=(nDArray-1)*nZeta + (nDArray-1)*nEta
      Call GetMem('MemR','ALLO','REAL',ipZeta,MemR)
      MemI=nZeta+nEta+2
      Call GetMem('MemI','ALLO','INTE',ipiZet,MemI)
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
      TCP3  = Zero
      rnint = Zero
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize memory for eval_int. Observe that semi-direct
*     calculations require the memory to be fix in between the
*     iterations. Hence, the memory pool must be handled externally in
*     those cases.
*
      If (XMem_Status.eq.Inactive) Call SetMem_Ints(0,0)
*
*     Call GetMem('S_I','Check','Real',iDum,iDum)
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
