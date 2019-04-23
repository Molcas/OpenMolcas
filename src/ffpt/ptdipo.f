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
      Subroutine PtDipo(H0,Ovlp,RR,nSize,Temp,nTemp)
*
************************************************************************
*                                                                      *
*     Objective: Construct the modified Hamiltonian,                   *
*                i.e., add dipole perturbation operator                *
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
      Data    Debug/.False./
      Dimension idum(1)
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Check if the command has been specified on input                 *
*                                                                      *
*----------------------------------------------------------------------*
*
      Call qEnter('PTDIPO')
*
      Exec=.false.
      Exec=Exec.or.ComStk(2,1,1,1)
      Exec=Exec.or.ComStk(2,1,1,2)
      Exec=Exec.or.ComStk(2,1,1,3)
      If ( .not.Exec ) then
         Call qExit ('PTDIPO')
         Return
      End If
*
*----------------------------------------------------------------------*
*     Loop over all components and add constributions to the           *
*     one-electron Hamiltonian                                         *
*----------------------------------------------------------------------*
*
      Do iComp=1,3
        If ( ComStk(2,1,1,iComp) ) Then
          Label='MltPl  1'
          iRc=-1
          iOpt1=1
          iOpt2=2
          iSyLbl=0
          Alpha=-ComVal(2,1,1,iComp)
          Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
          nInts=idum(1)
          If ( iRc.ne.0 ) Goto 991
          Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
          If ( iRc.ne.0 ) Goto 991
          Call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
          If ( Debug ) Then
            Write (6,'(6X,A,F8.6)') 'weight =',Alpha
            PriLbl='MltPl  1; Comp =    '
            Write(PriLbl(19:20),'(I2)') iComp
            Call PrDiOp(PriLbl,nSym,nBas,Temp)
            Write (6,*) 'Nuclear contribution=',Temp(nInts+4)
          End If
          call daxpy_(nInts,Alpha,Temp,1,H0,1)
          H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
        End If
      End Do
*
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
*
      Call qExit ('PTDIPO')
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
991   Write (6,*) 'PtDipi: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call QTrace
      Call Abend()
      End
