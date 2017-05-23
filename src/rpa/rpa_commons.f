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
      Subroutine RPA_SetInc()
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
      Integer i, j
      ! rpa_config
      Reference='Non'
      RPAModel='None@Non'
      DFTFunctional='Not defined     '
      dRPA=.false.
      SOSEX=.false.
      doCD=.false.
      doDF=.false.
      doLDF=.false.
      LumOrb=.false.
      iPrint=0
      ! rpa_data
      Do i=1,mTitle
         Write(Title(i),'(80A1)') (' ',j=1,80)
      End Do
      nTitle=0
      nSym=0
      Call iZero(nFreeze,2)
      Call iZero(nBas,8)
      Call iZero(nOrb,8)
      Call iZero(nFro,16)
      Call iZero(nDel,16)
      Call iZero(nOcc,16)
      Call iZero(nVir,16)
      Call iZero(ip_CMO,2)
      Call iZero(l_CMO,2)
      Call iZero(ip_EMO,2)
      Call iZero(l_EMO,2)
      Call iZero(ip_OccEn,2)
      Call iZero(l_OccEn,2)
      Call iZero(ip_VirEn,2)
      Call iZero(l_VirEn,2)
      NuclearRepulsionEnergy(1)=0.0d0
      End
************************************************************************
      Subroutine RPA_SetIntegralRepresentation()
      Implicit None
#include "rpa_config.fh"
      Call DecideOnCholesky(doCD)
      Call DecideOnDF(doDF)
      Call DecideOnLocalDF(doLDF)
      If (doLDF) Then
         doCD=.false.
         doDF=.false.
      Else If (doDF) Then
         doCD=.false.
         doLDF=.false.
      Else If (doCD) Then
         doDF=.false.
         doLDF=.false.
      End If
      End
************************************************************************
      Subroutine RPA_CheckIntegralRepresentation()
      Implicit None
#include "rpa_config.fh"
      If (.not.(doCD.or.doDF.or.doLDF)) Then
         Call RPA_Warn(2,'RPA requires CD, DF, or LDF. '//
     *                        'Conventional integrals not implemented.')
      End If
      End
************************************************************************
      Integer Function RPA_iUHF()
      Implicit None
#include "rpa_config.fh"
      Integer iUHF
      If (Reference(1:1).eq.'R') Then
         iUHF=1
      Else If (Reference(1:1).eq.'U') Then
         iUHF=2
      Else
         Write(6,'(A,A)') 'Reference=',Reference
         Call RPA_Warn(3,'Unable to determine iUHF in RPA')
         iUHF=-1
      End If
      RPA_iUHF=iUHF
      End
************************************************************************
      Integer Function RPA_LENIN4()
      Implicit None
#include "Molcas.fh"
      RPA_LENIN4=LENIN4
      End
