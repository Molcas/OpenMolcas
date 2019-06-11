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
      Subroutine PtGLbl(H0,Ovlp,RR,nSize,Temp,nTemp)
*
************************************************************************
*                                                                      *
*     Objective: Construct the modified Hamiltonian,                   *
*                i.e., add any perturbation defined by the GLBL input  *
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
      Logical Debug,Exec
      Data    Debug/.false./
      Dimension idum(1)
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Check if the command has been specified on input                 *
*                                                                      *
*----------------------------------------------------------------------*
*
      Call qEnter('PTGLBL')
*
      Exec=.false.
      Exec=Exec.or.ComStk(3,0,0,0)
      If ( .not.Exec ) then
         Call qExit ('PTGLBL')
         Return
      End If
*
*----------------------------------------------------------------------*
*     Loop over all components and add constributions to the           *
*     one-electron Hamiltonian                                         *
*----------------------------------------------------------------------*
*
      Do iLbl=1,mLbl
         Label=gLblN(iLbl)
         iComp=gLblC(iLbl)
         Alpha=gLblW(iLbl)
         If ( Debug ) Then
           Write(6,'(6X,5A,I2,2A,G9.2)')
     &     'GLBL    ',
     &     'label  ="',gLblN(iLbl),'",',
     &     'comp   =',gLblC(iLbl),',',
     &     'weight =',gLblW(iLbl)
         End If
         iRc=-1
         iOpt1=1
         iOpt2=2
         iSyLbl=0
         Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
         nInts=idum(1)
         If ( iRc.ne.0 ) Goto 991
         Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
         If ( iRc.ne.0 ) Goto 991
         Call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
         If ( Debug ) Then
           PriLbl=Label//'; Comp =    '
           Write(PriLbl(19:20),'(I2)') iComp
           Call PrDiOp(PriLbl,nSym,nBas,Temp)
         End If
         call daxpy_(nInts,Alpha,Temp,1,H0,1)
         H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
      End Do
*
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
*
      Call qExit ('PTGLBL')
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
991   Write (6,*) 'PtGlbl: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call QTrace
      Call Abend()
      End
