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
      Subroutine PtQuad(H0,Ovlp,RR,nSize,Temp,nTemp)
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
      Integer DiComp(3)
      Data    DiComp/1,4,6/
      Logical Debug,Exec,Orig,Diag
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
      Call qEnter('PtQuad')
*
      Exec=.false.
      Exec=Exec.or.ComStk(2,2,1,1)
      Exec=Exec.or.ComStk(2,2,1,2)
      Exec=Exec.or.ComStk(2,2,1,3)
      Exec=Exec.or.ComStk(2,2,1,4)
      Exec=Exec.or.ComStk(2,2,1,5)
      Exec=Exec.or.ComStk(2,2,1,6)
      Exec=Exec.or.ComStk(2,2,1,7)
      If ( .not.Exec ) then
         Call qExit('PtQuad')
         Return
      End If
*
*----------------------------------------------------------------------*
*     Check if a origin has been specified                             *
*     The unspecified components of the origin are set to 0.0 !        *
*     If no origin has been given pick the center of mass!             *
*----------------------------------------------------------------------*
*
c        Do i=1,nAtoms
c          Print *,i,(Coor(j,i),j=1,3)
c        End Do

      Orig=.false.
      Orig=Orig.or.ComStk(2,2,2,1)
      Orig=Orig.or.ComStk(2,2,2,2)
      Orig=Orig.or.ComStk(2,2,2,3)
      Orig=Orig.or.ComStk(2,2,2,4)
*
      If ( Orig ) Then
        XOrig=0.0
        YOrig=0.0
        ZOrig=0.0
        If ( ComStk(2,2,2,1) ) XOrig=ComVal(2,2,2,1)
        If ( ComStk(2,2,2,2) ) YOrig=ComVal(2,2,2,2)
        If ( ComStk(2,2,2,3) ) ZOrig=ComVal(2,2,2,3)
        If ( ComStk(2,2,2,4) ) Then
          iAtm=INT(ComVal(2,2,2,4))
          If ( Debug )
     *      Write(6,'(6X,A,I2)')'Origin of perturbation is centered '//
     *                          'at atom ',iAtm
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
*     Check if a diagonal component has been specified.                *
*     If so, compute R**2 first.                                       *
*----------------------------------------------------------------------*
*
      Diag=.false.
      Do iDiag=1,3
        Diag=Diag.or.ComStk(2,2,1,DiComp(iDiag))
      End Do
*-----or if RR option is used
      Diag=Diag.or.ComStk(2,2,1,7)
*
      If ( Diag ) Then
        Do iDiag=1,3
          Label='MltPl  2'
          iRc=-1
          iOpt1=1
          iOpt2=0
          iSyLbl=0
          iComp=DiComp(iDiag)
          Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
          nInts=idum(1)
          If ( iRc.ne.0 ) Goto 991
          Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
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
          If ( iComp.eq.1 ) Then
             call dcopy_(nInts,[0.0D0],0,RR,1)
             RR(nInts+4)=0.0D0
          End If
          Alpha=-0.5D0
          call daxpy_(nInts,Alpha,Temp,1,RR,1)
          RR(nInts+4)=RR(nInts+4)+Alpha*Temp(nInts+4)
          If ( Debug ) Then
            Write (6,'(6X,A,F8.6)') 'weight =',ComVal(2,2,1,iComp)
            PriLbl='MltPl  2; Comp =    '
            Write(PriLbl(19:20),'(I2)') iComp
            Call PrDiOp(PriLbl,nSym,nBas,Temp)
          End If
        End Do
        If ( Debug ) Then
          PriLbl='MltPl  2; Comp =R**2'
          Call PrDiOp(PriLbl,nSym,nBas,RR)
        End If
      End If
*
*----------------------------------------------------------------------*
*     Loop over components                                             *
*----------------------------------------------------------------------*
*
      jDiag=1
      iDiag=DiComp(jDiag)
      Do iComp=1,6
        If ( ComStk(2,2,1,iComp) ) Then
          Label='MltPl  2'
          iRc=-1
          iOpt1=1
          iOpt2=0
          iSyLbl=0
          Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
          nInts=idum(1)
          If ( iRc.ne.0 ) Goto 991
          Call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
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
          If ( iComp.eq.iDiag ) Then
            Alpha=-ComVal(2,2,1,iComp)
            call daxpy_(nInts,Alpha,RR,1,H0,1)
            H0(nInts+4)=H0(nInts+4)-Alpha*RR(nInts+4)
            Alpha=1.5*Alpha
            call daxpy_(nInts,Alpha,Temp,1,H0,1)
            H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
          Else
            Alpha=-1.5D0*ComVal(2,2,1,iComp)
            call daxpy_(nInts,Alpha,Temp,1,H0,1)
            H0(nInts+4)=H0(nInts+4)-Alpha*Temp(nInts+4)
            If ( Debug ) Then
              Write (6,'(6X,A,F8.6)') 'weight =',ComVal(2,2,1,iComp)
              PriLbl='MltPl  2; Comp =    '
              Write(PriLbl(19:20),'(I2)') iComp
              Call PrDiOp(PriLbl,nSym,nBas,Temp)
            End If
          End If
        End If
        If ( iComp.eq.iDiag .and. iComp.ne.DiComp(3)) Then
          jDiag=jDiag+1
          iDiag=DiComp(jDiag)
        End If
      End Do
*-----RR option
      If ( ComStk(2,2,1,7) ) Then
         Alpha=2.0D0*(ComVal(2,2,1,7))
         CALL DAXPY_(nInts,Alpha,RR,1,H0,1)
         H0(nInts+4)=H0(nInts+4)-Alpha*RR(nInts+4)
      End If
*
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
*
      Call qExit('PtQuad')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Ovlp)
*
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
*
991   Write (6,*) 'PtQuad: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call QTrace
      Call Abend()
      End
