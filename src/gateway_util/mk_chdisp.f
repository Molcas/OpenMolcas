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
* Copyright (C) 1991,1992,2008, Roland Lindh                           *
************************************************************************
      SubRoutine Mk_ChDisp()
************************************************************************
*                                                                      *
* Object: To generate the displacement labels
*                                                                      *
* Called from: Gateway
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September '91                                            *
*                                                                      *
*             Modified to complement GetInf, January '92.              *
************************************************************************
      use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
      Integer nDisp(0:7), DegDisp(MxAtom*3)
      Character ChDisp(MxAtom*3)*(LENIN6)
      Logical TstFnc
      Character*1 xyz(0:2)
      Data xyz/'x','y','z'/
*                                                                      *
************************************************************************
*                                                                      *
      LuWr=6
*
      nCnttp_Valence=0
      Do iCnttp = 1, nCnttp
         If (dbsc(iCnttp)%Aux) Exit
         nCnttp_Valence = nCnttp_Valence+1
      End Do
*
      mDisp = 0
      mdc = 0
      Do 10 iCnttp = 1, nCnttp_Valence
         If (dbsc(iCnttp)%pChrg) Then
             mdc = mdc + dbsc(iCnttp)%nCntr
             Go To 10
         End If
         Do iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
            mDisp = mDisp + 3*(nIrrep/dc(mdc)%nStab)
         End Do
 10   Continue
*                                                                      *
************************************************************************
*                                                                      *
      iDisp = 0
      Do iIrrep = 0, nIrrep-1
*        Loop over basis function definitions
         mdc = 0
         mc = 1
         nDisp(iIrrep)=0
         Do iCnttp = 1, nCnttp_Valence
*           Loop over unique centers associated with this basis set.
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
*              Loop over the cartesian components
               Do iCar = 0, 2
                  iComp = 2**iCar
                  If ( TstFnc(iOper,nIrrep,dc(mdc)%iCoSet,
     &                iChTbl,iIrrep,iComp,dc(mdc)%nStab) .and.
     &                .Not.dbsc(iCnttp)%pChrg ) Then
                     iDisp = iDisp + 1
                     ChDisp(iDisp)=' '
                     Write (ChDisp(iDisp)(1:(LENIN6)),'(A,1X,A1)')
     &                     dc(mdc)%LblCnt(1:LENIN4),xyz(iCar)
                     DegDisp(iDisp)=nIrrep/dc(mdc)%nStab
                     nDisp(iIrrep)=nDisp(iIrrep)+1
                  End If
*
               End Do
               mc = mc + nIrrep/dc(mdc)%nStab
            End Do
         End Do
*
      End Do
*
      If (iDisp.ne.mDisp) Then
         Call WarningMessage(2,
     &              ' Wrong number of symmetry adapted displacements')
         Write (6,*) iDisp,'=/=',mDisp
         Call Abend()
      End If
*
      Call Put_iScalar('nChDisp',iDisp)
      Call Put_cArray('ChDisp',ChDisp(1),(LENIN6)*iDisp)
      Call Put_iArray('nDisp',nDisp,nIrrep)
      Call Put_iArray('DegDisp',DegDisp,iDisp)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
