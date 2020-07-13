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
      SubRoutine BitMap_Localisation(PreFix)
      use Index_arrays, only: iSO2Sh
      Implicit Real*8 (a-h,o-z)
      Character*2 PreFix
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"

      Character*19 SecNam
      Parameter (SecNam = 'BitMap_Localisation')

      Integer iBas(8)
      Logical Indexation, DoF, DoG

C     Define iBas.
C     ------------

      nBasT   = nBas(1)
      iBas(1) = 0
      Do iSym = 2,nSym
         iBas(iSym) = nBasT
         nBasT = nBasT + nBas(iSym)
      End Do

C     Allocate and define some index arrays from Seward.
C     --------------------------------------------------

      DoF  = .false.
      nDiff = 0
      Call IniSew(DoF,nDiff)
      nShell = -1
      Indexation = .true.
      ThrAO = 0.0d0
      DoF = .false.
      DoG = .false.
      Call Setup_Ints(nShell,Indexation,ThrAO,DoF,DoG)
      If (nShell .lt. 1) Then
         Call SysAbendMsg(SecNam,'Setup_Ints failed!','nShell < 1')
      End If

C     Allocate max. sym. block of density matrix
C     and shell based density and CMO matrices.
C     ------------------------------------------

      MxBa = nBas(1)
      MxOr = nOrb2Loc(1)
      Do iSym = 2,nSym
         MxBa = max(MxBa,nBas(iSym))
         MxOr = max(MxOr,nOrb2Loc(iSym))
      End Do
      lDen = MxBa**2
      lDSh = nShell**2
      lCSh = nShell*MxOr
      lXSh = lCSh
      Call GetMem('BMpLoc','Allo','Real',ipDen,lDen)
      Call GetMem('Dsh','Allo','Real',ipDSh,lDSh)
      Call GetMem('Csh','Allo','Real',ipCSh,lCSh)
      Call GetMem('Xsh','Allo','Real',ipXSh,lXSh)

C     Compute density matrix, Den = CC^T, and set shell based matrices.
C     Generate bitmap and perform sparsity analysis.
C     -----------------------------------------------------------------

      kC = ipMOrig
      kX = ipCMO
      Do iSym = 1,nSym
         kC1 = kC + nBas(iSym)*nFro(iSym)
         Call GetDens_Localisation(Work(ipDen),Work(kC1),nBas(iSym),
     &                             nOrb2Loc(iSym))
         iOff = 1
         Call GetSh_Localisation(Work(ipDen),nBas(iSym),nBas(iSym),
     &                           Work(ipDSh),nShell,iSO2Sh(iOff),2,
     &                           AnaNrm)
         Call GetSh_Localisation(Work(kC1),nBas(iSym),nOrb2Loc(iSym),
     &                           Work(ipCSh),nShell,iSO2Sh(iOff),1,
     &                           AnaNrm)
         kX1 = kX + nBas(iSym)*nFro(iSym)
         Call GetSh_Localisation(Work(kX1),nBas(iSym),nOrb2Loc(iSym),
     &                           Work(ipXSh),nShell,iSO2Sh(iOff),1,
     &                           AnaNrm)
         Call GenBMp_Localisation(Work(ipDSh),Work(ipCSh),Work(ipXSh),
     &                            nShell,iSym,'r','r','r',PreFix)
         Call Anasize_Localisation(Work(ipDSh),Work(ipCSh),Work(ipXSh),
     &                             nShell,nOrb2Loc(iSym),iSym)
         n2 = nBas(iSym)**2
         kC = kC + n2
         kX = kX + n2
      End Do
      Write(6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

C     De-allocations.
C     ---------------

      Call GetMem('Xsh','Free','Real',ipXSh,lXSh)
      Call GetMem('Csh','Free','Real',ipCSh,lCSh)
      Call GetMem('Dsh','Free','Real',ipDSh,lDSh)
      Call GetMem('BMpLoc','Free','Real',ipDen,lDen)
      DoF = .false.
      DoG = .false.
      Call Term_Ints(DoF,DoG)

      End
