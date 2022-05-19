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
      Subroutine MMCount(natom,nAtMM,IsMM)
      Implicit Real*8 (A-H,O-Z)
*
*     Count the number of MM atoms
*
#include "espf.fh"
#include "stdalloc.fh"
      Integer, Intent(InOut):: nAtom
      Integer, Intent(InOut):: IsMM(nAtom)
      Integer, Intent(Out):: nAtMM
      Integer, Allocatable:: IsMM1(:), NTC(:)
*
      Logical Exist
*
      iPL = iPL_espf()
      ipIsMM = ip_iDummy
*
      Call Qpg_iArray('IsMM',Exist,nBla)
      If (.not.Exist) Then
         Write(6,'(A)') 'MMCount: IsMM not on the runfile'
         Call Abend()
      End If
      If (nBla.le.0) Then
         Write(6,'(A,I5)') 'MMCount: IsMM bad length:',nBla
         Call Abend()
      End If
      Call mma_allocate(IsMM1,nBla,Label='IsMM1')
      Call Get_iArray('IsMM',IsMM1,nBla)

      Call mma_allocate(NTC,natom,Label='NTC')
      Call Get_iArray('Atom -> Basis',NTC,natom)

      Do iAtom = 1, natom
         IsMM(iAtom) = IsMM1(NTC(iAtom))
      End Do

      Call mma_deallocate(NTC)
      Call mma_deallocate(IsMM1)
*
      nAtMM = 0
      Do iAt = 1, natom
         If (IsMM(iAt).eq.1) nAtMM = nAtMM + 1
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

      Return
      End
