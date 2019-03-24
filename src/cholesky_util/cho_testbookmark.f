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
      Subroutine Cho_TestBookmark(irc,verbose,is1CCD)
      Implicit None
      Integer irc
      Logical verbose
      Logical is1CCD
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character*16 SecNam
      Parameter (SecNam='Cho_TestBookmark')

      Logical dealloc

      Integer nVec(8)
      Integer jrc
      Integer iSym
      Integer NumChoBak(8)
      Integer NumChTBak
      Integer ipD
      Integer lD
      Integer ipDX
      Integer lDX

      Real*8  delta(8)
      Real*8  Thr
      Real*8  ThrComBak
      Real*8  ErrMx

      Logical Cho_1Center_Bak


      irc=0

      If (verbose) Then
         Call Cho_Head('Output from '//SecNam,'=',80,6)
      End If

      ! test 1: asking for accuracy below decomposition threshold
      !         should result in failure (return code 1)
      Thr=ThrCom*1.0d-1
      Call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
      If (jrc.eq.-1) Then ! bookmarks not available
         irc=-1
         If (verbose) Then
            Write(6,'(A,A)')
     &      'Cho_X_Bookmark returned -1 [not available].',
     &      'No further testing performed!'
         End If
         Return
      End If
      If (jrc.eq.1) Then
         If (verbose) Call Cho_TestBookmark_Prt(1,'passed')
      Else
         irc=irc+1
         If (verbose) Call Cho_TestBookmark_Prt(1,'failed')
      End If

      ! test 2: asking for negative accuracy
      !         should result in failure (return code 1)
      Thr=-1.0d-12
      Call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
      If (jrc.eq.1) Then
         If (verbose) Call Cho_TestBookmark_Prt(2,'passed')
      Else
         irc=irc+1
         If (verbose) Call Cho_TestBookmark_Prt(2,'failed')
      End If

      ! test 3: asking for decomposition threshold should
      !         give the total number of vectors.
      Thr=ThrCom
      Call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
      If (jrc.eq.0) Then
         Do iSym=1,nSym
            If (nVec(iSym).ne.NumCho(iSym)) Then
               jrc=jrc+1
            End If
            If (delta(iSym).gt.ThrCom) Then
               jrc=jrc+1
            End If
         End Do
         If (jrc.eq.0) Then
            If (verbose) Call Cho_TestBookmark_Prt(3,'passed')
         Else
            irc=irc+1
            If (verbose) Call Cho_TestBookmark_Prt(3,'failed')
         End If
      Else
         irc=irc+1
         If (verbose) Call Cho_TestBookmark_Prt(3,'failed')
      End If

      ! test 4: a threshold above the decomposition threshold
      !         should result in fewer vectors and less accuracy.
      If (ThrCom.lt.1.0d-4) Then
         Thr=max(ThrCom*1.0d3,1.0d-14)
      Else
         Thr=max(ThrCom*1.0d2,1.0d-14)
      End If
      Call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
      If (jrc.eq.0) Then
         Do iSym=1,nSym
            If (delta(iSym).gt.Thr) Then
               jrc=jrc+1
            End If
            If (nVec(iSym).gt.NumCho(iSym)) Then
               jrc=jrc+1
            Else If (nVec(iSym).eq.NumCho(iSym)) Then
               If (delta(iSym).gt.ThrCom) Then
                  jrc=jrc+1
               End If
            End If
         End Do
         If (jrc.eq.0) Then
            If (verbose) Call Cho_TestBookmark_Prt(4,'passed')
         Else
            irc=irc+1
            If (verbose) Call Cho_TestBookmark_Prt(4,'failed')
         End If
      Else
         irc=irc+1
         If (verbose) Call Cho_TestBookmark_Prt(4,'failed')
      End If

      ! Test 5: check diagonal (only if previous tests passed).
      If (irc.ne.0) Then
         If (verbose) Call Cho_TestBookmark_Prt(5,'not executed')
      Else
         Cho_1Center_Bak=Cho_1Center
         Cho_1Center=is1CCD
         Do iSym=1,nSym
            NumChoBak(iSym)=NumCho(iSym)
         End Do
         NumChTBak=NumChT
         ThrComBak=ThrCom
         NumChT=0
         ThrCom=0.0d0
         Do iSym=1,nSym
            NumCho(iSym)=nVec(iSym)
            NumChT=NumChT+NumCho(iSym)
            ThrCom=max(ThrCom,delta(iSym))
         End Do
         lD=nnBstRT(1)
         Call GetMem('BkmDia','Allo','Real',ipD,lD)
         Call Cho_X_CalcChoDiag(jrc,Work(ipD))
         If (jrc.eq.0) Then
            lDX=nnBstRT(1)
            Call GetMem('BkmXDia','Allo','Real',ipDX,lDX)
            Call Cho_IODiag(Work(ipDX),2)
            Call dAXPY_(nnBstRT(1),-1.0d0,Work(ipD),1,Work(ipDX),1)
            If (Cho_1Center) Then
               Call Cho_TestBookmark_1CInit(dealloc)
               Call Cho_MaxAbsDiag(Work(ipDX),1,ErrMx)
               If (dealloc) Call Cho_TestBookmark_1CFinal()
            Else
               Call Cho_MaxAbsDiag(Work(ipDX),1,ErrMx)
            End If
            Call FZero(DiaMax,nSym)
            Call FZero(DiaMaxT,nSym)
            Call GetMem('BkmXDia','Free','Real',ipDX,lDX)
            If (abs(ErrMx-ThrCom).gt.1.0d-12) Then
               irc=irc+1
               If (verbose) Call Cho_TestBookmark_Prt(5,'failed')
            Else
               If (verbose) Call Cho_TestBookmark_Prt(5,'passed')
            End If
         Else
            irc=irc+1
            If (verbose) Call Cho_TestBookmark_Prt(5,'failed')
         End If
         Call GetMem('BkmDia','Free','Real',ipD,lD)
         Do iSym=1,nSym
            NumCho(iSym)=NumChoBak(iSym)
         End Do
         NumChT=NumChTBak
         ThrCom=ThrComBak
         Cho_1Center=Cho_1Center_Bak
      End If

      End
      Subroutine Cho_TestBookmark_Prt(TestNumber,PassFail)
      Implicit None
      Integer TestNumber
      Character*(*) PassFail
      Write(6,'(A,I3,1X,A)') 'Test',TestNumber,PassFail
      End
      Subroutine Cho_TestBookmark_1CInit(AllocatedHere)
      Implicit None
      Logical AllocatedHere
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer irc

      If (l_iAtomShl.lt.1) Then
         l_iAtomShl=nShell
         Call GetMem('iAtomShl','Allo','Inte',ip_iAtomShl,l_iAtomShl)
         irc=-1
         Call Cho_SetAtomShl(irc,iWork(ip_iAtomShl),l_iAtomShl)
         If (irc.ne.0) Then
            Write(6,'(A,I4)')
     &      'Cho_TestBookmark_1Cinit: Cho_SetAtomShl returned',irc
            Call Cho_Quit('shell-to-atom init failed!',104)
         End If
         AllocatedHere=.True.
      Else
         AllocatedHere=.False.
      End If

      End
      Subroutine Cho_TestBookmark_1CFinal()
      Implicit None
#include "choptr.fh"

      If (l_iAtomShl.gt.0) Then
         Call GetMem('iAtomShl','Free','Inte',ip_iAtomShl,l_iAtomShl)
         ip_iAtomShl=0
         l_iAtomShl=0
      End If

      End
