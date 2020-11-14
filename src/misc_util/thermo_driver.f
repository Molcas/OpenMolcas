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
      Subroutine Thermo_Driver(UserT,UserP,nUserPT,nsRot,
     &                              EVal,nFreq,lSlapaf)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Integer nUserPT, nsRot, nFreq
      Real*8 UserT(64), UserP, Eval(*)
      Logical lSlapaf ! If .True. then Thermo_Driver called by SLAPAF
      Logical lTest
      lTest = .False.
*
      If (lSlapaf) Then
         Call Get_iScalar('NSYM',nSym)
         If (nSym.ne.1) Then
            Write (6,'(A)') 'WARNING: '
     &                //'No thermochemistry analysis conducted for '
     &           //'numerical frequencies unless no symmetry is used!'
            Return
         End If
      End If
*
      Write(6,*)
      Call CollapseOutput(1,'Thermochemistry')
      Write(6,*)
      Write(6,'(1X,A)')'*********************'
      Write(6,'(1X,A)')'*                   *'
      Write(6,'(1X,A)')'*  THERMOCHEMISTRY  *'
      Write(6,'(1X,A)')'*                   *'
      Write(6,'(1X,A)')'*********************'
      Write(6,*)
*
      If (lTest) then
        Write(6,*)'----------------------------------------------------'
        Write(6,*)'[Thermo_Driver] Input Data:'
        Write(6,*)'    UserP=',UserP,'  nsRot=',nsRot,'nUserPT=',nUserPT
        Write(6,*)'    UserT(1-5)==',(UserT(i),i=1,5)
        Write(6,'(A,I3,A,256F8.2)')'  nFreq=',nFreq,
     &            '  Freq(i)==',(EVal(i),i=1,nFreq)
        Write(6,*)'----------------------------------------------------'
        Call XFlush(6)
      EndIf
*
      Call Rotation(TotalM,TRotA,TRotB,TRotC,nsRot,nAtom,lSlapaf)
      Call Get_iScalar('Multiplicity',iMult)
*
      If (lTest) then
        Write(6,*)' Calling ThermoChem,  iMult=',iMult
        Write(6,*)' UserP=',UserP,'  nsRot=',nsRot,'  nAtom=',nAtom
        Write(6,*)' TotalM,TRotA,TRotB,TRotC==',TotalM,TRotA,TRotB,TRotC
        Write(6,'(A,I3,A,256F8.2)')' nFreq=',nFreq,
     &            '  Freq(i)==',(EVal(i),i=1,nFreq)
        Call XFlush(6)
      EndIf
*
      Call ThermoChem_(UserT,UserP,TotalM,TRotA,TRotB,TRotC,nUserPT,
     &                 nsRot,iMult,nAtom,EVal,nFreq,lSlapaf)
      Call CollapseOutput(0,'Thermochemistry')

      Return
      End
