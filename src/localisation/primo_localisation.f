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
*               2011, Francesco Aquilante                              *
************************************************************************
      Subroutine PriMO_Localisation(Header,PrOcc,PrEne,ThrOcc,ThrEne,
     &                              nSym,nBas,nOrb,Nme,Ene,Occ,CMO,
     &                              iPrForm,IndxT)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Print MOs.
C     This is a work-around to fix bugs when orbitals are deleted.
C
C     F. Aquilante, Nov 2011   (Print only non-deleted orbitals)
C
      Implicit None
#include "Molcas.fh"
      Character*(*) Header
      Character*(LENIN8) Nme(*)
      Logical PrOcc, PrEne
      Real*8  ThrOcc, ThrEne
      Integer nSym
      Integer nBas(nSym), nOrb(nSym)
      Real*8  Ene(*), Occ(*), CMO(*)
      Integer iPrForm
      Integer IndxT(*)
#include "WrkSpc.fh"

      Integer ip_CMO, l_CMO
      Integer ip_Occ, l_OCc
      Integer ip_EOr, l_EOr
      Integer nOrbT, iSym, k1, k2
      Integer k, kk, ik, nOrb_(8)

      Call Icopy(nSym,nOrb,1,nOrb_,1)
      kk=0
      Do iSym=1,nSym
         Do k=1,nBas(iSym)
            ik=kk+k
            If (IndxT(ik).eq.7) nOrb(iSym)=nOrb(iSym)-1
         End Do
         kk=kk+nBas(iSym)
      End Do

      nOrbT=nOrb(1)
      Do iSym=2,nSym
         nOrbT=nOrbT+nOrb(iSym)
      End Do

      l_CMO=nBas(1)*nOrb(1)
      Do iSym=2,nSym
         l_CMO=l_CMO+nBas(iSym)*nOrb(iSym)
      End Do
      l_Occ=nOrbT
      l_EOr=nOrbT

      Call GetMem('CMO_','Allo','Real',ip_CMO,l_CMO)
      Call GetMem('Occ_','Allo','Real',ip_Occ,l_Occ)
      Call GetMem('Eor_','Allo','Real',ip_EOr,l_EOr)

      k1=1
      k2=ip_CMO
      Do iSym=1,nSym
         Call dCopy_(nBas(iSym)*nOrb(iSym),CMO(k1),1,Work(k2),1)
         k1=k1+nBas(iSym)*nBas(iSym)
         k2=k2+nBas(iSym)*nOrb(iSym)
      End Do

      k1=1
      k2=ip_Occ
      Do iSym=1,nSym
         Call dCopy_(nOrb(iSym),Occ(k1),1,Work(k2),1)
         k1=k1+nBas(iSym)
         k2=k2+nOrb(iSym)
      End Do

      k1=1
      k2=ip_EOr
      Do iSym=1,nSym
         Call dCopy_(nOrb(iSym),Ene(k1),1,Work(k2),1)
         k1=k1+nBas(iSym)
         k2=k2+nOrb(iSym)
      End Do

      Call PriMO(Header,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb,
     &           Nme,Work(ip_EOr),Work(ip_Occ),Work(ip_CMO),iPrForm)

      Call GetMem('Eor_','Free','Real',ip_EOr,l_EOr)
      Call GetMem('Occ_','Free','Real',ip_Occ,l_Occ)
      Call GetMem('CMO_','Free','Real',ip_CMO,l_CMO)

      Call Icopy(nSym,nOrb_,1,nOrb,1)

      End
