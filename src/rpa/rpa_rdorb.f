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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RPA_RdOrb()
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Read orbitals and orbital energies from InpOrb or from Runfile.
C
      Implicit None
#include "rpa_config.fh"

      Character*9 SecNam
      Parameter (SecNam='RPA_RdOrb')


      If (LumOrb) Then
         ! read from InpOrb
         Call RPA_RdOrb_FromInpOrb()
      Else
         ! read from Runfile
         Call RPA_RdOrb_FromRunfile()
      End If


      End
************************************************************************
      Subroutine RPA_RdOrb_FromInpOrb()
      Implicit None

      Character*20 SecNam
      Parameter (SecNam='RPA_RdOrb_FromInpOrb')

      Call RPA_Warn(3,
     *    SecNam//': Reading orbitals from INPORB not implemented yet')

      End
************************************************************************
      Subroutine RPA_RdOrb_FromRunfile()
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Character*21 SecNam
      Parameter (SecNam='RPA_RdOrb_FromRunfile')

      Integer  RPA_iUHF
      External RPA_iUHF

      Integer iUHF
      Integer iSym
      Integer i
      Integer nB, nB2
      Integer ip, ipO, ipV

      ! Restricted (1) or unrestricted (2)
      iUHF=RPA_iUHF()

      ! Allocate memory for CMO
      l_CMO(1)=nBas(1)*nOrb(1)
      nB2=nBas(1)**2
      Do iSym=2,nSym
         l_CMO(1)=l_CMO(1)+nBas(iSym)*nOrb(iSym)
         nB2=nB2+nBas(iSym)**2
      End Do
      Call GetMem('CMO(RPA)','Allo','Real',ip_CMO(1),l_CMO(1))
      If (iUHF.eq.2) Then
         l_CMO(2)=l_CMO(1)
         Call GetMem('CMO(RPA)','Allo','Real',ip_CMO(2),l_CMO(2))
      Else
         ip_CMO(2)=0
         l_CMO(2)=0
      End If

      ! Read CMO array(s) from Runfile
      Call Get_CMO(Work(ip_CMO(1)),nB2)
      If (iUHF.eq.2) Then
         Call Get_dArray('CMO_ab',Work(ip_CMO(2)),nB2)
      End If

      ! Allocate memory for orbital energies
      nB=nBas(1)
      Do iSym=2,nSym
         nB=nB+nBas(iSym)
      End Do
      Do i=1,iUHF
         l_OccEn(i)=nOcc(1,i)
         l_VirEn(i)=nVir(1,i)
         Do iSym=2,nSym
            l_OccEn(i)=l_OccEn(i)+nOcc(iSym,i)
            l_VirEn(i)=l_VirEn(i)+nVir(iSym,i)
         End Do
         Call GetMem('OccEn','Allo','Real',ip_OccEn(i),l_OccEn(i))
         Call GetMem('VirEn','Allo','Real',ip_VirEn(i),l_VirEn(i))
      End Do
      If (iUHF.eq.1) Then
         ip_OccEn(2)=0
         l_OccEn(2)=0
         ip_VirEn(2)=0
         l_VirEn(2)=0
      End If

      ! Read orbital energies from Runfile
      Call Get_OrbE(ip_EMO(1),l_EMO(1))
      If (l_EMO(1).ne.nB) Then
         Call RPA_Warn(3,SecNam//': unexpected EMO dimension')
      End If
      ip=ip_EMO(1)
      ipO=ip_OccEn(1)
      ipV=ip_VirEn(1)
      Do iSym=1,nSym
         Call dCopy_(nOcc(iSym,1),Work(ip),1,Work(ipO),1)
         Call dCopy_(nVir(iSym,1),Work(ip+nOcc(iSym,1)),1,Work(ipV),1)
         ip=ip+nOrb(iSym)
         ipO=ipO+nOcc(iSym,1)
         ipV=ipV+nVir(iSym,1)
      End Do
      If (iUHF.eq.2) Then
         l_EMO(2)=l_EMO(1)
         Call GetMem('EMO(RPA)','Allo','Real',ip_EMO(2),l_EMO(2))
         Call Get_dArray('OrbE_ab',Work(ip_EMO(2)),l_EMO(2))
         ip=ip_EMO(2)
         ipO=ip_OccEn(2)
         ipV=ip_VirEn(2)
         Do iSym=1,nSym
            Call dCopy_(nOcc(iSym,2),Work(ip),1,Work(ipO),1)
           Call dCopy_(nVir(iSym,2),Work(ip+nOcc(iSym,2)),1,Work(ipV),1)
            ip=ip+nOrb(iSym)
            ipO=ipO+nOcc(iSym,2)
            ipV=ipV+nVir(iSym,2)
         End Do
      End If

      End
