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
      Subroutine PtOkt0(H0,Ovlp,RR,nSize,Temp,nTemp)
*
************************************************************************
*                                                                      *
*     Objective: Construct the modified Hamiltonian,                   *
*                i.e., add quadrupole perturbation operator            *
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
      Dimension Cntr(3)
      Logical Debug,Exec,Orig
      Data    Debug/.False./
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Check if the command has been specified on input                 *
*                                                                      *
*----------------------------------------------------------------------*
*
      Call qEnter('PtOkt0')
*
      Exec=.false.
      Exec=Exec.or.ComStk(2,6,1,1)
      Exec=Exec.or.ComStk(2,6,1,2)
      Exec=Exec.or.ComStk(2,6,1,3)
      Exec=Exec.or.ComStk(2,6,1,4)
      Exec=Exec.or.ComStk(2,6,1,5)
      Exec=Exec.or.ComStk(2,6,1,6)
      Exec=Exec.or.ComStk(2,6,1,7)
      Exec=Exec.or.ComStk(2,6,1,8)
      Exec=Exec.or.ComStk(2,6,1,9)
      Exec=Exec.or.ComStk(2,6,1,10)
      If ( .not.Exec ) then
         Call qExit('PtOkt0')
         Return
      End If
*
*----------------------------------------------------------------------*
*     Check if a origin has been specified                             *
*     The unspecified components of the origin are set to 0.0 !        *
*     If no origin has been given pick the center of mass!             *
*----------------------------------------------------------------------*
*
      Orig=.false.
      Orig=Orig.or.ComStk(2,6,2,1)
      Orig=Orig.or.ComStk(2,6,2,2)
      Orig=Orig.or.ComStk(2,6,2,3)
      Orig=Orig.or.ComStk(2,6,2,4)
*
      If ( Orig ) Then
        XOrig=0.0
        YOrig=0.0
        ZOrig=0.0
        If ( ComStk(2,6,2,1) ) XOrig=ComVal(2,6,2,1)
        If ( ComStk(2,6,2,2) ) YOrig=ComVal(2,6,2,2)
        If ( ComStk(2,6,2,3) ) ZOrig=ComVal(2,6,2,3)
        If ( ComStk(2,6,2,4) ) Then
          iAtm=INT(ComVal(2,6,2,4))
          If ( iAtm.lt.0 .or. iAtm.gt.nAtoms ) Then
             Write (6,*) 'PtOkt0: You specified a invalid atom number'
     &                 //' as the origin of the perturbation operator.'
             Call QTrace
             Call Abend()
          End If
          XOrig=Coor(1,iAtm)
          YOrig=Coor(2,iAtm)
          ZOrig=Coor(3,iAtm)
        End If
      Else
        Call Get_dArray('Center of Mass',Cntr,3)
        XOrig=Cntr(1)
        YOrig=Cntr(2)
        ZOrig=Cntr(3)
      End If
      If ( Debug )
     *  Write(6,'(6X,A,3F12.6)')'Origin of the perturbation operator =',
     *  XOrig,YOrig,ZOrig
*
*----------------------------------------------------------------------*
*     Loop over components                                             *
*----------------------------------------------------------------------*
*
      Do iComp=1,10
        If ( ComStk(2,6,1,iComp) ) Then
          If ( iComp.eq.1 .or. iComp.eq.4 .or. iComp.eq.6 ) Then
            Call PtOkt1('XR',Temp,RR)
          Else If ( iComp.eq.2 .or. iComp.eq.7 .or. iComp.eq.9 ) Then
            Call PtOkt1('YR',Temp,RR)
          Else If ( iComp.eq.3 .or. iComp.eq.8 .or. iComp.eq.10 ) Then
            Call PtOkt1('ZR',Temp,RR)
          End If
          Label='MltPl  3'
          PriLbl='MltPl  3; Comp =    '
          Write(PriLbl(19:20),'(I2)') iComp
          iRc=-1
          iOpt1=1
          iOpt2=0
          iSyLbl=0
          Call iRdOne(iRc,iOpt1,Label,iComp,nInts,iSyLbl)
          If ( iRc.ne.0 ) Goto 991
          Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
          If ( iRc.ne.0 ) Goto 991
          Call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
          X=Temp(nInts+1)
          Y=Temp(nInts+2)
          Z=Temp(nInts+3)
          If ( X.ne.XOrig .or. Y.ne.YOrig .or. Z.ne.ZOrig ) Then
             Write (6,*) 'PtOkt0: Input error, no matching center'
     &                 //' is found.'
             Call QTrace
             Call Abend()
          End If
          Alpha=5.0d0
          Call DSCAL_(nInts+4,Alpha,Temp,1)
          If ( iComp.eq.1 .or. iComp.eq.7 .or. iComp.eq.10 ) Then
            Alpha=-3.0d0
            call daxpy_(nInts,Alpha,RR,1,Temp,1)
            Temp(nInts+4)=Temp(nInts+4)+Alpha*RR(nInts+4)
          Else If ( iComp.ne.5 ) Then
            Alpha=-1.0d0
            call daxpy_(nInts,Alpha,RR,1,Temp,1)
            Temp(nInts+4)=Temp(nInts+4)+Alpha*RR(nInts+4)
          End If
          Alpha=0.5d0*ComVal(2,6,1,iComp)
          call daxpy_(nInts,Alpha,Temp,1,H0,1)
          H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
          If ( Debug ) Then
            Write (6,'(6X,A,F8.6)') 'weight =',ComVal(2,6,1,iComp)
            PriLbl='MltPl  3; Comp =    '
            Write(PriLbl(19:20),'(I2)') iComp
            Call PrDiOp(PriLbl,nSym,nBas,Temp)
          End If
        End If
      End Do
*
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
*
      Call qExit('PtOkt0')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Ovlp)
*
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
*
991   Write (6,*) 'PtOkt0: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call QTrace
      Call Abend()
      End
