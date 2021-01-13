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
      Subroutine PrintQ(rK,qLbl,nq,nQQ,LuIC,rMult)
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
      Real*8 rK(nq,nQQ), rMult(nq)
      Character*14  qLbl(nq), Line*80, filnam*16
      Logical Start
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      iPrint=99
#else
      iRout=122
      iPrint=nPrint(iRout)
#endif
*
      Lu=6
      Thr=0.001D+00 ! Threshold for printout.
*
      If (iPrint.le.5) Go To 99
      Write (Lu,*)
      Call CollapseOutput(1,'Internal coordinate section')
*
      If (iPrint.ge.6) Then
         Write (Lu,*)
         Write (Lu,'(80A)') ('*',i=1,80)
         Write (Lu,*) ' Auto-Defined Internal coordinates'
         Write (Lu,'(80A)') ('-',i=1,80)
         Write (Lu,'(A)') '  Primitive Internal Coordinates:'
      Else
         Write (Lu,*)
         Write (Lu,'(A)') '  Redundant Internal Coordinates:'
         Write (Lu,*)
      End If
*
      Rewind(LuIC)
 999  Continue
         Read (LuIC,'(A)',END=998) Line
         Write (Lu,'(A)') Line
         Go To 999
 998  Continue
      Rewind(LuIC)
      If (iPrint.lt.6) Go To 99
*
      Write (Lu,'(A)') '  Internal Coordinates:'
      Do iQQ = 1, nQQ
         Write(Line,'(A,I3.3,A)') 'q',iQQ,' ='
         iF=7
         jq=0
         Start=.True.
         Do iq = 1, nq
            temp=Abs(rK(iq,iQQ))
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
                  Write(Line(iF:iE),'(A,F10.8,4A)') ' ',rK(iq,iQQ),
     &                                            ' ',qLbl(iq)(1:4),' '
               Else
                  iE=iF+17
                  Write(Line(iF:iE),'(A,F10.8,4A)') '+ ',rK(iq,iQQ),
     &                                            ' ',qLbl(iq)(1:4),' '
               End If
               iF=iE+1
            End If
         End Do
         Write (Lu,'(A)') Line
      End Do
      Write (Lu,'(80A)') ('*',i=1,80)
      Call CollapseOutput(0,'Internal coordinate section')
 99   Continue
*
*     Write linear combinations to disc
*
      LuTmp=11
      filnam='SPCINX'
      call molcas_binaryopen_vanilla(luTmp,filnam)
c     Open(luTmp,File=filnam,Form='unformatted',Status='unknown')
      ReWind (LuTmp)
*
*---- put in degeneracy factor so that ddot will work.
*
      Write (LuTmp) nq,nQQ
      Do iq = 1, nq
         Write (LuTmp) qLbl(iq),(rMult(iq)*rK(iq,iQQ),iQQ=1,nQQ)
      End Do
*
      Close  (LuTmp)
*
      If (iPrint.ge.10.and.nQQ.le.12) Then
         Write (Lu,*)
         Write (Lu,*) ' Nonredundant internal coordinates'
      End If
      If (iPrint.ge.6.and.nQQ.le.12) Then
         Write (Lu,*)
         Write (Lu,*) ' Number of redundant coordinates:',nq
         Write (Lu,*)
      End If
      If (iPrint.ge.10.and.nQQ.le.12) Then
         Write (Lu,'(A,E10.3)') ' Threshold for printout:',Thr
         IncQQ = 8
         Do iiQQ = 1, nQQ, IncQQ
            mQQ=Min(nQQ,iiQQ+IncQQ-1)
            Write (Lu,*)
            Write (Lu,'(14X,8I10)') (iQQ,iQQ=iiQQ,mQQ)
            Do iq = 1, nq
               temp=Sqrt(DDot_(nQQ,rK(iq,1),nq,rK(iq,1),nq))
               If (temp.gt.Thr)
     &            Write (Lu,'(A,8F10.6)')
     &                  qLbl(iq),(rK(iq,iQQ),iQQ=iiQQ,mQQ)
            End Do
            Write (Lu,*)
         End Do
      End If
*
      Return
      End
