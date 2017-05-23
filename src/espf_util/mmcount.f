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
      Subroutine MMCount(natom,nAtMM,ipIsMM)
      Implicit Real*8 (A-H,O-Z)
*
*     Count the number of MM atoms
*
#include "espf.fh"
*
      Logical Exist
*
      Call QEnter('MMCount')
      iPL = iPL_espf()
      ipIsMM = ip_iDummy
*
      Call Qpg_iArray('IsMM',Exist,nBla)
      If (.not.Exist) Then
         Write(6,'(A)') 'MMCount: IsMM not on the runfile'
         Call Abend
      End If
      If (nBla.le.0) Then
         Write(6,'(A,I5)') 'MMCount: IsMM bad length:',nBla
         Call Abend
      End If
      Call GetMem('Is MM','Allo','Inte',ipIsMM1,nBla)
      Call Get_iArray('IsMM',iWork(ipIsMM1),nBla)
      Call GetMem('AtoToBas','Allo','Inte',ipNTC,natom)
      Call Get_iArray('Atom -> Basis',iWork(ipNTC),natom)
      Call GetMem('IsMM for atoms','Allo','Inte',ipIsMM,natom)
      Do iAtom = 0, natom-1
         iWork(ipIsMM+iAtom) = iWork(ipIsMM1+iWork(ipNTC+iAtom)-1)
      End Do
      Call GetMem('AtoToBas','Free','Inte',ipNTC,natom)
      Call GetMem('Is MM','Free','Inte',ipIsMM1,nBla)
*
      nAtMM = 0
      Do iAt = 1, natom
         If (iWork(ipIsMM+iAt-1).eq.1) nAtMM = nAtMM + 1
      End Do
      If (nAtMM.lt.0) Then
         Write(6,'(A)') 'Error in MMCount: nAtMM < 0!'
         Call Quit_OnUserError()
      Else If (nAtMM.gt.natom) Then
         Write(6,'(A)') 'Error in MMCount: nAtMM >= natom!'
         Call Quit_OnUserError()
      Else If(nAtMM.ne.0.and.iPL.ge.3) Then
         Write(6,'(A,I5,A)') ' QM/MM: found ',nAtMM,' MM atoms'
      End If
      Call QExit('MMCount')
      Return
      End
