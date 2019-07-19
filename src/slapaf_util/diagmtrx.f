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
      Subroutine DiagMtrx(H,nH,iNeg)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      character*16 filnam
      Real*8 H(nH,nH)
      Logical Exist
*

*
      Call QEnter('DiagMtrx')
      Lu=6
      iRout=21
      iPrint=nPrint(iRout)
*
      Call GetMem('EVal','Allo','Real',ipEVal,nH*(nH+1)/2)
      Call GetMem('EVec','Allo','Real',ipEVec,nH*nH)
*
*---- Copy elements for H
*
      SumHii=Zero
      Do i = 1, nH
         Do j = 1, i
            ij=i*(i-1)/2 + j + ipEval -1
            Work(ij)=H(i,j)
         End Do
         SumHii=SumHii+H(i,i)
      End Do
*     Write (Lu,*) ' SumHii=',SumHii
*
*---- Set up a unit matrix
*
      call dcopy_(nH*nH,[Zero],0,Work(ipEVec),1)
      call dcopy_(nH,[One],0,Work(ipEVec),nH+1)
*
*---- Compute eigenvalues and eigenvectors
*
      Call NIDiag_new(Work(ipEVal),Work(ipEVec),nH,nH,0)
      Call Jacord(Work(ipEVal),Work(ipEVec),nH,nH)
*
*---- Print out the result
*
      iNeg=0
      Do i = 1, nH
         If (Work(i*(i+1)/2+ipEVal-1).lt.Zero) iNeg=iNeg+1
      End Do
      IF (iprint.gt.5) THEN
        Write (Lu,*)
        Write (Lu,*)'************************************************'//
     &              '*****************'
        Write (Lu,*)'* Eigenvalues and Eigenvectors of the Hessian   '//
     &              '                *'
        Write (Lu,*)'************************************************'//
     &              '*****************'
      END IF
*
      filnam='SPCINX'
      call f_Inquire(filnam,Exist)
*
      If (Exist .AND. iprint.gt.5) Then
*
*        Read linear combinations from disc
*

         LuTmp=11
         call molcas_binaryopen_Vanilla(luTmp,filnam)
c         Open(luTmp,File=filnam,Form='unformatted',Status='unknown')
         ReWind (LuTmp)
*
         Read (LuTmp) nq,nQQ
*
         If (nQQ.ne.nH) then
           Write (Lu,*)
           Write (Lu,*) 'Eigenvalues of the Hessian'
           Write (Lu,*)
           Write (Lu,'(1X,10F10.5)') (Work(i*(i+1)/2+ipEVal-1),i=1,nH)
         EndIf
         If (nQQ.eq.nH) Then
*
           Call GetMem('rK','Allo','Real',iprK,nq*nQQ)
           Call GetMem('qEVec','Allo','Real',ipqEVec,nq*nH)
*
           Call Print_qEVec(Work(ipEVec),nH,ipEVal,nq,
     &                      Work(iprK),Work(ipqEVec),LuTmp)
*
           Call GetMem('qEVec','Free','Real',ipqEVec,nq*nH)
           Call GetMem('rK','Free','Real',iprK,nq*nQQ)
*
         Else
*
           Write (Lu,*)
           Write (Lu,*) 'Eigenvectors of the Hessian'
           Write (Lu,*)
           Do i = 1, nH
              Write (Lu,'(1X,10F10.5)')
     &              (Work((j-1)*nH+i+ipEVec-1),j=1,nH)
           End Do
         End If
*
         Close (LuTmp)
*
      Else If (iprint.gt.5) Then
*
        Call Print_qEVec2(nH,ipEVal,Work(ipEVec))
*
      End If
*
      Call GetMem('EVec','Free','Real',ipEVec,nH*nH)
      Call GetMem('EVal','Free','Real',ipEVal,nH*(nH+1)/2)
*
      Call QExit('DiagMtrx')
      Return
      End


      Subroutine Print_qEVec(EVec,nH,ipEVal,nq,rK,qEVec,LuTmp)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 EVec(nH,nH), rK(nq,nH), qEVec(nq,nH)
      Character*14 qLbl(nq)
*
      Do iq = 1, nq
         Read (LuTmp) qLbl(iq),(rK(iq,iQQ),iQQ=1,nH)
      End Do
*
      Call DGEMM_('N','N',
     &            nq,nH,nH,
     &            1.0d0,rK,nq,
     &            EVec,nH,
     &            0.0d0,qEVec,nq)
*
      Lu=6
      Thr=0.0001D0
      IncQQ = 5
      Do iiQQ = 1, nH, IncQQ
         mQQ=Min(nH,iiQQ+IncQQ-1)
         Write (Lu,*)
         Write (Lu,'(14X,5I10)') (iQQ,iQQ=iiQQ,mQQ)
         Write (Lu,'(1X,A,5F10.6)') 'Eigenvalues   ',
     &               (Work(iQQ*(iQQ+1)/2+ipEVal-1),iQQ=iiQQ,mQQ)
         Write (Lu,*)
         Do iq = 1, nq
            temp=Sqrt(DDot_(nH,qEVec(iq,1),nq,qEVec(iq,1),nq))
            If (temp.gt.Thr)
     &         Write (Lu,'(1X,A,5F10.6)')
     &               qLbl(iq),(qEVec(iq,iQQ),iQQ=iiQQ,mQQ)
         End Do
         Write (Lu,*)
      End Do
*
      Return
      End

      Subroutine Print_qEVec2(nH,ipEVal,EVec)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 EVec(nH,nH)
      Character*14 qLbl(nH)
      Character*14 cLbl
      Character*120 Temp
*
* --- Skip Primitive Coords
*
      Lu_UDIC=91
      Temp='UDIC'
      Call molcas_open(Lu_UDIC,Temp)
  10  Read(Lu_UDIC,'(A)') Temp
      Call UpCase(Temp)
      If (Temp(1:4).EQ.'VARY') Go To 20
      GoTo 10
*
* --- Read Internal Coords Labels
*
   20 Do iLines = 1, nH
   40    Read(Lu_UDIC,'(A)') Temp
         Call UpCase(Temp)
         If (Temp(1:3).EQ.'FIX') Go To 40
         cLbl = ' '
         Do j = 1, 14
           If(Temp(j:j).EQ.' '.or.Temp(j:j).EQ.'=') GoTo 30
           cLbl(j:j)=Temp(j:j)
         EndDo
  30     Continue
         qLbl(iLines) = cLbl
      EndDo

      Lu=6
      IncQQ = 5
      Do iiQQ = 1, nH, IncQQ
         mQQ=Min(nH,iiQQ+IncQQ-1)
         Write (Lu,*)
         Write (Lu,'(14X,5I10)') (iQQ,iQQ=iiQQ,mQQ)
         Write (Lu,'(1X,A,5F10.6)') 'Eigenvalues   ',
     &               (Work(iQQ*(iQQ+1)/2+ipEVal-1),iQQ=iiQQ,mQQ)
         Write (Lu,*)
         Do iq = 1, nH
             Write (Lu,'(1X,A,5F10.6)') qLbl(iq),
     &               (EVec(iq,iQQ),iQQ=iiQQ,mQQ)
         End Do
         Write (Lu,*)
      End Do
*
      Close(Lu_UDIC)
      Return
      End
