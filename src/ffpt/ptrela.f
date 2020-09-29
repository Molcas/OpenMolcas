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
      Subroutine PtRela(H0,Ovlp,RR,nSize,Temp,nTemp)
*
************************************************************************
*                                                                      *
*     Objective: Add the relativitic perturbation operator to          *
*                the one-electron Hamiltonian                          *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*

#include "input.fh"
*
      Real*8 H0(nSize), Ovlp(nSize), RR(nSize), Temp(nTemp)
      Character*8 Label
      Character*20 PriLbl
      Logical Debug
      Data    Debug/.false./
      Dimension idum(1)
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Load mass-velocity and one-electron Darwin contact term.         *
*     Add them to the one-electron Hamiltonian.                        *
*                                                                      *
*----------------------------------------------------------------------*
*
*
      If ( .not.ComStk(2,5,0,1) ) then
         Return
      End If
*
      Label='MassVel '
      iRc=-1
      iOpt1=1
      iOpt2=2
      iSyLbl=0
      iComp=1
      Alpha=ComVal(2,5,0,1)
      Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
      nInts=idum(1)
      If ( iRc.ne.0 ) Goto 991
      Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
      If ( iRc.ne.0 ) Goto 991
      call daxpy_(nInts,Alpha,Temp,1,H0,1)
      H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
      If ( Debug ) Then
        Write (6,'(6X,A,F8.6)') 'weight =',Alpha
        PriLbl='Mass-Velocity term  '
        Write(PriLbl(19:20),'(I2)') iComp
        Call PrDiOp(PriLbl,nSym,nBas,Temp)
      End If
      Label='Darwin  '
      iRc=-1
      iOpt1=1
      iOpt2=2
      iSyLbl=0
      iComp=1
      Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
      nInts=idum(1)
      If ( iRc.ne.0 ) Goto 991
      Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
      If ( iRc.ne.0 ) Goto 991
      call daxpy_(nInts,Alpha,Temp,1,H0,1)
      H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
      If ( Debug ) Then
        Write (6,'(6X,A,F8.6)') 'weight =',Alpha
        PriLbl='1el. Darwin term    '
        Write(PriLbl(19:20),'(I2)') iComp
        Call PrDiOp(PriLbl,nSym,nBas,Temp)
      End If
*
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Ovlp)
        Call Unused_real_array(RR)
      End If
*
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
*
991   Write (6,*) 'PtRela: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call QTrace
      Call Abend()
      End
