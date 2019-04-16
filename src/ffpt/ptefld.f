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
      Subroutine PtEfld(H0,Ovlp,RR,nSize,Temp,nTemp)
*
************************************************************************
*                                                                      *
*     Objective: Construct the modified Hamiltonian                    *
*                i.e., add electric field perturbation operator        *
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
      Logical Debug,Exec,Orig,NoCntr
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
      Call qEnter('PtEfld')
*
      Exec=.false.
      Exec=Exec.or.ComStk(2,3,1,1)
      Exec=Exec.or.ComStk(2,3,1,2)
      Exec=Exec.or.ComStk(2,3,1,3)
      If ( .not.Exec ) then
         Call qExit('PtEfld')
         Return
      End If
*
*----------------------------------------------------------------------*
*     Check if a origin has been specified                             *
*     The unspecified components of the origin are set to 0.0 !        *
*----------------------------------------------------------------------*
*
      Orig=.false.
      Orig=Orig.or.ComStk(2,3,2,1)
      Orig=Orig.or.ComStk(2,3,2,2)
      Orig=Orig.or.ComStk(2,3,2,3)
      Orig=Orig.or.ComStk(2,3,2,4)
      If ( .not.Orig ) Then
         Write (6,*) 'PtElfd: No matching center is found.'
         Call QTrace
         Call Abend()
      End If
*
      XOrig=0.0
      YOrig=0.0
      ZOrig=0.0
      If ( ComStk(2,3,2,1) ) XOrig=ComVal(2,3,2,1)
      If ( ComStk(2,3,2,2) ) YOrig=ComVal(2,3,2,2)
      If ( ComStk(2,3,2,3) ) ZOrig=ComVal(2,3,2,3)
      If ( ComStk(2,3,2,4) ) Then
        iAtm=INT(ComVal(2,3,2,4))
        If ( iAtm.lt.0 .or. iAtm.gt.nAtoms ) Then
           Write (6,*) 'PtEfld: You specified a invalid atom number as'
     &               //' the origin of the perturbation operator.'
           Call QTrace
           Call Abend()
        End If
        XOrig=Coor(1,iAtm)
        YOrig=Coor(2,iAtm)
        ZOrig=Coor(3,iAtm)
      End If
      If ( Debug )
     *  Write(6,'(6X,A,3F12.6)')'Origin of perturbation operator =',
     *  XOrig,YOrig,ZOrig
*
*----------------------------------------------------------------------*
*     Loop over the max possible number of centers and                 *
*     search for coincidence in the origin definitions                 *
*----------------------------------------------------------------------*
*
      MxCntr=9999
      NoCntr=.true.
      Do iCntr=1,MxCntr
        If ( NoCntr ) Then
          Label='EF1     '
          Write(Label(4:8),'(I5)')iCntr
          iRc=-1
          iOpt1=1
          iOpt2=2
          iSyLbl=0
          Do iComp=1,3
            Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
            nInts=idum(1)
            If ( iRc.eq.0 ) Then
              Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
              X=Temp(nInts+1)
              Y=Temp(nInts+2)
              Z=Temp(nInts+3)
              If ( X.eq.XOrig .and. Y.eq.YOrig .and. Z.eq.ZOrig )
     *          NoCntr=.false.
            End If
          End Do
        End If
      End Do
      If ( NoCntr ) Then
         Write (6,*) 'PtEfld: You missed to specify the origin of '
     &             //'the operator.'
         Call QTrace
         Call Abend()
      End If
      If ( Debug )
     *  Write(6,'(6X,A,A)')'Label of perturbation operator =', Label
*
*----------------------------------------------------------------------*
*     If centers match read the integrals and accumulate contribution  *
*----------------------------------------------------------------------*
*
      Do iComp=1,3
        If ( ComStk(2,3,1,iComp) ) Then
          iRc=-1
          iOpt1=1
          iOpt2=2
          iSyLbl=0
          Alpha=-ComVal(2,3,1,iComp)
          Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
          nInts=idum(1)
          If ( iRc.ne.0 ) Goto 991
          Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
          If ( iRc.ne.0 ) Goto 991
          Call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
          call daxpy_(nInts,Alpha,Temp,1,H0,1)
          H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
          If ( Debug ) Then
            Write (6,'(6X,A,F8.6)') 'weight =',Alpha
            PriLbl='        ; Comp =    '
            PriLbl(1:8)=Label
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
      Call qExit('PtEfld')
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
991   Write (6,*) 'PtEfld: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call QTrace
      Call Abend()
      End
