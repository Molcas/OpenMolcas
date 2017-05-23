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
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      SubRoutine List2(Title,Lbl,gq,nAtom,nInter,Smmtrc)
************************************************************************
*                                                                      *
* Object: to print cartesian internal coordinates.                     *
*                                                                      *
* Called from: RlxCtl                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             1993                                                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
*
      Real*8 gq(3*nAtom,nInter)
      Character Lbl(nAtom)*(*), Title*(*)
      Logical Smmtrc(3*nAtom)
*
      n_qLbl=3*nAtom
      nChar=4*n_qLbl
      Call GetMem('qLbl','Allo','Char',ip_qLbl,nChar)
      Call List2_(Title,Lbl,gq,nAtom,nInter,Smmtrc,cWork(ip_qLbl),
     &            n_qLbl)
      Call GetMem('qLbl','Free','Char',ip_qLbl,nChar)
*
      Return
      End
      SubRoutine List2_(Title,Lbl,gq,nAtom,nInter,Smmtrc,qLbl,n_qLbl)
************************************************************************
*                                                                      *
* Object: to print cartesian internal coordinates.                     *
*                                                                      *
* Called from: RlxCtl                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             1993                                                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
*
      Real*8 gq(3*nAtom,nInter)
      Character Lbl(nAtom)*(*), Format*72, Title*(*), Line*80
      Character*4 qLbl(n_qLbl), qLbl_tmp*14
      Logical Start, Smmtrc(3*nAtom)
      character*16 filnam
      Lu=6
*
      iRout = 119
      iPrint = nPrint(iRout)
      Call qEnter('List')
*
      Thr=0.001D+00 ! Threshold for printout.
*
      mInt=3*nAtom
*
      Write (Lu,*)
      Call CollapseOutput(1,'Internal coordinates')
      Write (Lu,*)
      Write (Lu,*) ' Specification of the internal coordinates '
     &          //'according to the user-defined internal'
      Write(Lu,*) ' coordinate format.'
      Write (Lu,*)
      Write (Lu,'(A)') 'Internal Coordinates'
      iq = 0
      Do igq = 1, mInt, 3
         If (Smmtrc(igq  )) Then
            iq=iq+1
            Write (qLbl(igq  ),'(A,I3.3)') 'c',iq
            Write (Lu,'(3A)')
     &             qLbl(igq  ),' = Cartesian x ',Lbl((igq+2)/3)
         End If
         If (Smmtrc(igq+1)) Then
            iq=iq+1
            Write (qLbl(igq+1),'(A,I3.3)') 'c',iq
            Write (Lu,'(3A)')
     &             qLbl(igq+1),' = Cartesian y ',Lbl((igq+2)/3)
         End If
         If (Smmtrc(igq+2)) Then
            iq=iq+1
            Write (qLbl(igq+2),'(A,I3.3)') 'c',iq
            Write (Lu,'(3A)')
     &             qLbl(igq+2),' = Cartesian z ',Lbl((igq+2)/3)
         End If
      End Do
      Write (Lu,'(A)') 'Vary'
      Do iQQ = 1, nInter
         Write(Line,'(A,I3.3,A)') 'q',iQQ,' ='
         iF=7
         jq=0
         Start=.True.
         Do iq = 1, mInt
            temp=Abs(gq(iq,iQQ))
            If (temp.gt.Thr) Then
               jq = jq + 1
               If (jq.gt.4) Then
                  Line(80:80)='&'
                  Write (Lu,'(A)') Line
                  Line=' '
                  iF=6
                  jq = 1
                  Start=.False.
               End If
               If (jq.eq.1.and.Start) Then
                  iE=iF+16
                  Write(Line(iF:iE),'(A,F10.8,4A)') ' ',gq(iq,iQQ),
     &                                            ' ',qLbl(iq),' '
               Else
                  iE=iF+17
                  Write(Line(iF:iE),'(A,F10.8,4A)') '+ ',gq(iq,iQQ),
     &                                            ' ',qLbl(iq),' '
               End If
               iF=iE+1
            End If
         End Do
         Write (Lu,'(A)') Line
      End Do
      Write (Lu,'(A)') 'End Of Internal Coordinates'
      Call CollapseOutput(0,'Internal coordinates')
*
*     Write linear combinations to disc
*
      LuTmp=11
      filnam='SPCINX'
      call molcas_binaryopen_vanilla(luTmp,filnam)
c      Open(luTmp,File=filnam,Form='unformatted',Status='unknown')
      ReWind (LuTmp)
*
      Write (LuTmp) mInt,nInter
      Do iq = 1, mInt
         qLbl_tmp=qLbl(iq)
         Write (LuTmp) qLbl_tmp,(gq(iq,iQQ),iQQ=1,nInter)
      End Do
*
      Close  (LuTmp)
*
      Write (Lu,*)
      Call CollapseOutput(1,Title)
*
      MxWdth=132
      nLbl=8+1
      nRow=9
      inc = Min((MxWdth-nLbl)/nRow,nInter)
*
      Do 10 ii = 1, nInter, inc
         Write (Lu,*)
         Write(Format,'(A,I2,A)') '(A,1X,',inc,'(I5,4X))'
         Write (Lu,Format) 'Internal',(i,i=ii,Min(ii+inc-1,nInter))
         Write (Lu,*)
         Write(Format,'(A,I2,A)') '(A4,A4,1X,',inc,'(F8.5,1X))'
         Do 20 igq = 1, mInt, 3
            Write (Lu,Format) Lbl((igq+2)/3),' x  ',
     &            (gq(igq  ,i),i=ii,Min(ii+inc-1,nInter))
            Write (Lu,Format) Lbl((igq+2)/3),' y  ',
     &            (gq(igq+1,i),i=ii,Min(ii+inc-1,nInter))
            Write (Lu,Format) Lbl((igq+2)/3),' z  ',
     &            (gq(igq+2,i),i=ii,Min(ii+inc-1,nInter))
 20      Continue
         Write (Lu,*)
 10   Continue
      Call CollapseOutput(0,Title)
*
      Call qExit('List')
      Return
      End
