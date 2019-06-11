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
      Subroutine PtOkt1(Oper,Temp1,Temp2)
*
************************************************************************
*                                                                      *
*     Objective: Construct the perturbation operator of the form       *
*                <X*R**2>=<XXX>+<XYY>+<XZZ>                            *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*

#include "input.fh"
*
      Dimension Temp1(*),Temp2(*)
*
      Character*2 Oper
      Character*8 Label
      Character*20 PriLbl
      Dimension Cntr(3)
      Integer xrComp(3)
      Data    xrComp/1,4,6/
      Integer yrComp(3)
      Data    yrComp/2,7,9/
      Integer zrComp(3)
      Data    zrComp/3,8,10/
      Logical Debug,Orig
      Data    Debug/.False./
      Dimension idum(1)
*
*----------------------------------------------------------------------*
*     Start procedure                                                  *
*     Check if a origin has been specified                             *
*     The unspecified components of the origin are set to 0.0 !        *
*     If no origin has been given pick the center of mass!             *
*----------------------------------------------------------------------*
*
      Call qEnter('PtOkt1')
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
      Do iComp=1,3
         If ( Oper.eq.'XR' ) jComp=xrComp(iComp)
         If ( Oper.eq.'YR' ) jComp=yrComp(iComp)
         If ( Oper.eq.'ZR' ) jComp=zrComp(iComp)
         Label='MltPl  3'
         PriLbl='MltPl  3; Comp =    '
         Write(PriLbl(19:20),'(I2)') jComp
         iRc=-1
         iOpt1=1
         iOpt2=0
         iSyLbl=0
         Call iRdOne(iRc,iOpt1,Label,jComp,idum,iSyLbl)
         nInts=idum(1)
         If ( iRc.ne.0 ) Goto 991
         Call RdOne(iRc,iOpt2,Label,jComp,Temp1,iSyLbl)
         Call CmpInt(Temp1,nInts,nBas,nSym,iSyLbl)
         X=Temp1(nInts+1)
         Y=Temp1(nInts+2)
         Z=Temp1(nInts+3)
         If ( X.ne.XOrig .or. Y.ne.YOrig .or. Z.ne.ZOrig ) Then
             Write (6,*) 'PtOkt1: Input error, no matching center'
     &                 //' is found.'
             Call QTrace
             Call Abend()
         End If
         If ( iComp.eq.1 ) Then
            call dcopy_(nInts,Temp1,1,Temp2,1)
            Temp2(nInts+4)=Temp1(nInts+4)
         Else
            Alpha=1.0d0
            call daxpy_(nInts,Alpha,Temp1,1,Temp2,1)
            Temp2(nInts+4)=Temp2(nInts+4)+Alpha*Temp1(nInts+4)
         End If
      End Do
*
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
*
      Call qExit('PtOkt1')
      Return
*
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
*
991   Write (6,*) 'PtOkt1: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call QTrace
      Call Abend()
      End
