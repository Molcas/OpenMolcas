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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RdVec_Localisation(nSym,nBas,nOrb,IndT,CMO,Occ,EOrb,
     &                              FName)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Read orbital info and return in a format suitable for module
C     localisation: deleted orbitals are included (as zeros). This is a
C     work-around to fix bugs when orbitals are deleted.
C
      Implicit None
      Integer nSym ! number of irreps
      Integer nBas(nSym)  ! number of basis functions
      Integer nOrb(nSym)  ! number of orbitals
      Integer IndT(*)  ! type indices, dim: nBas
      Real*8  CMO(*)   ! MO coefficients, dim: nBas*nBas
      Real*8  Occ(*)   ! Occupation numbers, dim: nBas
      Real*8  EOrb(*)  ! EOrb, dim: nBas
      Character*(*) FName ! filename with input orbitals
#include "warnings.fh"
#include "WrkSpc.fh"

      Character*18 SecNam
      Parameter (SecNam='RdVec_Localisation')

      Character*80 VTitle

      Integer nBasT, nOrbT, iSym
      Integer ip_CMO, l_CMO
      Integer ip_Occ, l_Occ
      Integer ip_EOr, l_EOr
      Integer ip_Ind, l_Ind
      Integer Lu, iUHF, iWarn, iErr, iWFType, k1, k2, i
      Integer mylen
      External mylen

      Real*8 Dummy(1)

      nBasT=nBas(1)
      nOrbT=nOrb(1)
      Do iSym=2,nSym
         nBasT=nBasT+nBas(iSym)
         nOrbT=nOrbT+nOrb(iSym)
      End Do

      l_CMO=nBas(1)*nOrb(1)
      Do iSym=2,nSym
         l_CMO=l_CMO+nBas(iSym)*nOrb(iSym)
      End Do
      l_Occ=nOrbT
      l_EOr=nOrbT
      l_Ind=nBasT
      Call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
      Call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
      Call GetMem('EOr_','Allo','Real',ip_EOr,l_EOr)
      Call GetMem('Ind_','Allo','Inte',ip_Ind,l_Ind)

      Lu=75
      iUHF=0  ! restricted HF
      iWarn=2 ! abend if nBas/nOrb info is inconsistent
      iErr=-1 ! init return code
      iWFType=-1 ! init wave function type
      Dummy(1)=9.9d9 ! dummy variable
      Call RdVec_(FName,Lu,'COEI',iUHF,nSym,nBas,nOrb,
     &            Work(ip_CMO),Dummy,Work(ip_Occ),Dummy,
     &            Work(ip_EOr),Dummy,iWork(ip_Ind),VTitle,iWarn,iErr,
     &            iWFType)
      If (iErr.ne.0) Then
         Call WarningMessage(2,
     &                     SecNam//': Non-zero return code from RdVec_')
         Write(6,'(A,A,I9)') SecNam,': RdVec_ returned code',iErr
         Call xFlush(6)
         Call xQuit(_RC_IO_ERROR_READ_)
      End If
      Write (6,*)
      Write (6,'(A)') ' Header from vector file:'
      Write (6,*)
      Write (6,'(A)') VTitle(:mylen(VTitle))
      Write (6,*)

      k1=ip_CMO
      k2=1
      Do iSym=1,nSym
         Call dCopy_(nBas(iSym)*nOrb(iSym),Work(k1),1,CMO(k2),1)
         Call Cho_dZero(CMO(k2+nBas(iSym)*nOrb(iSym)),
     &                  nBas(iSym)*(nBas(iSym)-nOrb(iSym)))
         k1=k1+nBas(iSym)*nOrb(iSym)
         k2=k2+nBas(iSym)*nBas(iSym)
      End Do

      k1=ip_Occ
      k2=1
      Do iSym=1,nSym
         Call dCopy_(nOrb(iSym),Work(k1),1,Occ(k2),1)
         Call Cho_dZero(Occ(k2+nOrb(iSym)),
     &                  nBas(iSym)-nOrb(iSym))
         k1=k1+nOrb(iSym)
         k2=k2+nBas(iSym)
      End Do

      k1=ip_EOr
      k2=1
      Dummy(1)=9.9d9
      Do iSym=1,nSym
         Call dCopy_(nOrb(iSym),Work(k1),1,EOrb(k2),1)
         Call dCopy_(nBas(iSym)-nOrb(iSym),Dummy(1),0,
     &                                    EOrb(k2+nOrb(iSym)),1)
         k1=k1+nOrb(iSym)
         k2=k2+nBas(iSym)
      End Do

      k1=ip_Ind
      k2=1
      Do iSym=1,nSym
         Call iCopy(nOrb(iSym),iWork(k1),1,IndT(k2),1)
         Do i=nOrb(iSym),nBas(iSym)-1
            IndT(k2+i)=7
         End Do
         k1=k1+nOrb(iSym)
         k2=k2+nBas(iSym)
      End Do

      Call GetMem('Ind_','Free','Inte',ip_Ind,l_Ind)
      Call GetMem('EOr_','Free','Real',ip_EOr,l_EOr)
      Call GetMem('Occ_','Free','Real',ip_Occ,l_Occ)
      Call GetMem('CMO_','Free','Real',ip_CMO,l_CMO)

      End
