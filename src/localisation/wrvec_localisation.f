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
      Subroutine WrVec_Localisation(FName,Lu,Label,nSym,nBas,nOrb,CMO,
     &                              Occ,EOrb,IndT,Title)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Write orbital info.
C     This is a work-around to fix bugs when orbitals are deleted.
C
      Implicit None
      Character*6  FName
      Integer Lu
      Character*(*) Label
      Integer nSym
      Integer nBas(nSym)
      Integer nOrb(nSym)
      Real*8  CMO(*)
      Real*8  Occ(*)
      Real*8  EOrb(*)
      Integer IndT(*)
      Character*(*) Title
#include "WrkSpc.fh"

      Integer ip_CMO, l_CMO
      Integer ip_Occ, l_Occ
      Integer ip_EOr, l_EOr
      Integer ip_Ind, l_Ind
      Integer iSym, k1, k2

      Logical Write_CMO, Write_Occ, Write_EOr, Write_Ind

      Write_CMO=index(Label,'C').ne.0
      Write_Occ=index(Label,'O').ne.0
      Write_EOr=index(Label,'E').ne.0
      Write_Ind=index(Label,'I').ne.0

      If (Write_CMO) Then
         l_CMO=nBas(1)*nOrb(1)
         Do iSym=2,nSym
            l_CMO=l_CMO+nBas(iSym)*nOrb(iSym)
         End Do
         Call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
         k1=1
         k2=ip_CMO
         Do iSym=1,nSym
            Call dCopy_(nBas(iSym)*nOrb(iSym),CMO(k1),1,Work(k2),1)
            k1=k1+nBas(iSym)*nBas(iSym)
            k2=k2+nBas(iSym)*nOrb(iSym)
         End Do
      Else
         l_CMO=1
         Call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
         Work(ip_CMO)=0.0d0
      End If

      If (Write_Occ) Then
         l_Occ=nOrb(1)
         Do iSym=2,nSym
            l_Occ=l_Occ+nOrb(iSym)
         End Do
         Call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
         k1=1
         k2=ip_Occ
         Do iSym=1,nSym
            Call dCopy_(nOrb(iSym),Occ(k1),1,Work(k2),1)
            k1=k1+nBas(iSym)
            k2=k2+nOrb(iSym)
         End Do
      Else
         l_Occ=1
         Call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
         Work(ip_Occ)=0.0d0
      End If

      If (Write_EOr) Then
         l_EOr=nOrb(1)
         Do iSym=2,nSym
            l_EOr=l_EOr+nOrb(iSym)
         End Do
         Call GetMem('EOr_','Allo','Real',ip_EOr,l_EOr)
         k1=1
         k2=ip_EOr
         Do iSym=1,nSym
            Call dCopy_(nOrb(iSym),EOrb(k1),1,Work(k2),1)
            k1=k1+nBas(iSym)
            k2=k2+nOrb(iSym)
         End Do
      Else
         l_EOr=1
         Call GetMem('EOr_','Allo','Real',ip_EOr,l_EOr)
         Work(ip_EOr)=0.0d0
      End If

      If (Write_Ind) Then
         l_Ind=56
         Call GetMem('Ind_','Allo','Inte',ip_Ind,l_Ind)
         Call iCopy(l_Ind,IndT,1,iWork(ip_Ind),1)
      Else
         l_Ind=1
         Call GetMem('Ind_','Allo','Inte',ip_Ind,l_Ind)
         iWork(ip_Ind)=0
      End If

      Call WrVec(FName,Lu,Label,nSym,nBas,nOrb,Work(ip_CMO),
     &           Work(ip_Occ),Work(ip_EOr),iWork(ip_Ind),Title)

      Call GetMem('Ind_','Free','Inte',ip_Ind,l_Ind)
      Call GetMem('EOr_','Free','Real',ip_EOr,l_EOr)
      Call GetMem('Occ_','Free','Real',ip_Occ,l_Occ)
      Call GetMem('CMO_','Free','Real',ip_CMO,l_CMO)

      End
