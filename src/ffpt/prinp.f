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
      Subroutine PrInp
*
************************************************************************
*                                                                      *
*     Objective: Write the title page on the standard output unit      *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*

#include "input.fh"
*
      Character*120 PrLine,BlLine,StLine
      Character*72  Data,Line
      Character*8 Fmt1,Fmt2
      Character*4 Com,Sub1,Sub2,Parm
      Integer StrnLn
      Logical clear, nice
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*                                                                      *
*----------------------------------------------------------------------*
*
      Call qEnter('PRINP')
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(PrLine)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      lPaper=132
      left=(lPaper-lLine)/2
      Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*
*     Print the project title                                          *
*----------------------------------------------------------------------*
      If ( mTit.gt.0 ) then
         Write(6,*)
         nLine=mTit+5
         Do i=1,nLine
           PrLine=BlLine
           If ( i.eq.1 .or. i.eq.nLine )
     &     PrLine=StLine
           If ( i.eq.3 )
     &     PrLine='Project:'
           If ( i.ge.4 .and. i.le.nLine-2 ) then
              PrLine=Title(i-3)
           End If
           Call Center(PrLine)
           Write(6,Fmt1) '*'//PrLine//'*'
         End Do
         Write(6,*)
      End If
*
*----------------------------------------------------------------------*
*     Print file identifier                                            *
*----------------------------------------------------------------------*
*
      Write(6,*)
      Write(6,'(6X,A)') 'Header of the ONEINT file:'
      Write(6,'(6X,A)') '--------------------------'
      Write(6,*)
      Write(Line,'(72A1)') (Header(i),i=1,72)
      Write(6,'(6X,A)') Line(:StrnLn(Line))
      Write(Line,'(72A1)') (Header(i),i=73,144)
      Write(6,'(6X,A)') Line(:StrnLn(Line))
      Write(6,*)
*
*----------------------------------------------------------------------*
*     Print the coordinates of the system                              *
*----------------------------------------------------------------------*
*
      Call PrCoor
*
*----------------------------------------------------------------------*
*     Print comand list                                                *
*----------------------------------------------------------------------*
*
      Write(6,*)
      Write(6,'(6X,A)')'The following tasks will be executed:'
      Write(6,'(6X,A)')'-------------------------------------'
      Write(6,*)
      Do iPr=1,120
         PrLine(iPr:iPr)=' '
      End Do
      Com='FFPT'
      nSub1=ComCtl(2,0,0)
      clear=.false.
      Do iSub1=1,nSub1
         Sub1=ComTab(2,iSub1,0,0)
         nSub2=ComCtl(2,iSub1,0)
         nParm=0
         If ( nSub2.eq.0 ) nParm=ComCtl(2,iSub1,1)
         ist=25
         Do iParm=0,nParm
            If ( ComStk(2,iSub1,0,iParm) ) Then
               Parm=ComTab(2,iSub1,0,iParm)
               PrLine(1:4)=Com
               PrLine(9:12)=Sub1
               PrLine(ist:ist+3)=Parm
               ist=ist+4
               z=ComVal(2,iSub1,0,iParm)
               Write(PrLine(ist:ist+7),'(F8.6)')z
               ist=ist+10
            End If
         End Do
         If ( PrLine(1:4).ne.'    ' ) Then
            If ( clear ) Then
               Write(6,'(14X,A)') PrLine(9:120)
            Else
               Write(6,'(6X,A)') PrLine
            End If
            Do iPr=1,120
               PrLine(iPr:iPr)=' '
            End Do
            clear=.true.
         End If
         clear=.false.
         Do iSub2=1,nSub2
            Sub2=ComTab(2,iSub1,iSub2,0)
            nParm=ComCtl(2,iSub1,iSub2)
            ist=25
            Do iParm=0,nParm
               If ( ComStk(2,iSub1,iSub2,iParm) ) Then
                  Parm=ComTab(2,iSub1,iSub2,iParm)
                  PrLine(1:4)=Com
                  PrLine(9:12)=Sub1
                  PrLine(17:20)=Sub2
                  PrLine(ist:ist+3)=Parm
                  ist=ist+4
                  z=ComVal(2,iSub1,iSub2,iParm)
                  Write(Data,'(F10.6)')z
                  Call LeftAd(Data)
                  iPoint=Index(Data,'.')
                  nice=.true.
                  Do i=iPoint+1,Strnln(Data)
                     If ( Data(i:i).ne.'0' ) nice=.false.
                  End Do
                  If ( nice ) then
                     PrLine(ist:ist+iPoint-1)=Data(1:1+iPoint-2)
                     ist=ist+iPoint+1
                  Else
                     PrLine(ist:ist+iPoint+5)=Data(1:1+iPoint+5)
                     ist=ist+iPoint+7
                  End If
               End If
            End Do
            If ( PrLine(1:4).ne.'    ' ) Then
               If ( clear ) Then
                  Write(6,'(22X,A)') PrLine(17:120)
               Else
                  Write(6,'(6X,A)') PrLine
               End If
               Do iPr=1,120
                  PrLine(iPr:iPr)=' '
               End Do
               clear=.true.
            End If
         End Do
      End Do
      Do iTbl=1,mLbl
         Write(6,'(6X,5A,I2,2A,F9.6)')
     &   'GLBL    ',
     &   'label="',gLblN(iTbl),'",',
     &   'comp=',gLblC(iTbl),',',
     &   'weight=',gLblW(iTbl)
      End Do
      Write(6,*)
*
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*
      Call qEXit('PRINP')
      Return
      End
