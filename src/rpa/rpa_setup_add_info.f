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
      Subroutine RPA_Setup_Add_Info()
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Data for RPA tests, checking the setup.
C     Testing these data is as much a test of the preceding SCF run as
C     a test of the RPA setup.
C
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Character*18 SecNam
      Parameter (SecNam='RPA_Setup_Add_Info')

      Integer  RPA_iUHF, Cho_X_GetTol
      External RPA_iUHF, Cho_X_GetTol
      Real*8   Cho_dSumElm, dDot_
      External Cho_dSumElm, dDot_

      Character*13 orbitals

      Integer Tol
      Integer iUHF
      Integer l_orbitals
      Integer iSpin
      Integer iSym
      Integer ipO, ipV
      Integer i

      Real*8 Tst(8)

      ! Check that molecular geometry is the expected one:
      ! nuclear repulsion energy
      Tol=12
      Call Add_Info('PotNuc',NuclearRepulsionEnergy,1,Tol)

      ! Check orbital spaces: sums and norms of orbital energies.
      Tol=min(Cho_X_GetTol(2),2)
      iUHF=RPA_iUHF()
      If (iUHF.eq.1) Then
         orbitals=' orbital'
         l_orbitals=8
      Else If (iUHF.eq.2) Then
         orbitals=' spin-orbital'
         l_orbitals=13
      Else
         Write(6,'(A,I6)') 'iUHF=',iUHF
         Call RPA_Warn(3,SecNam//': iUHF error')
         orbitals=' '
         l_orbitals=1
      End If
      Call fZero(Tst,8)
      Do iSpin=1,iUHF
         ipO=ip_OccEn(iUHF)
         ipV=ip_VirEn(iUHF)
         Do iSym=1,nSym
            Tst(1)=Tst(1)+Cho_dSumElm(Work(ipO),nFro(iSym,iSpin))
            Tst(2)=Tst(2)+dDot_(nFro(iSym,iSpin),Work(ipO),1,
     *                                          Work(ipO),1)
            ipO=ipO+nFro(iSym,iSpin)
            Tst(3)=Tst(3)+Cho_dSumElm(Work(ipO),nOcc(iSym,iSpin))
            Tst(4)=Tst(4)+dDot_(nOcc(iSym,iSpin),Work(ipO),1,
     *                                          Work(ipO),1)
            ipO=ipO+nOcc(iSym,iSpin)
            Tst(5)=Tst(5)+Cho_dSumElm(Work(ipV),nVir(iSym,iSpin))
            Tst(6)=Tst(6)+dDot_(nVir(iSym,iSpin),Work(ipV),1,
     *                                          Work(ipV),1)
            ipV=ipV+nVir(iSym,iSpin)
            Tst(7)=Tst(7)+Cho_dSumElm(Work(ipV),nDel(iSym,iSpin))
            Tst(8)=Tst(8)+dDot_(nDel(iSym,iSpin),Work(ipV),1,
     *                                          Work(ipV),1)
            ipV=ipV+nDel(iSym,iSpin)
         End Do
      End Do
      Do i=2,8,2
         Tst(i)=sqrt(Tst(i))
      End Do
      Call Add_Info(Reference//orbitals(1:l_orbitals)//' energy',Tst,8,
     *              Tol)

      End
