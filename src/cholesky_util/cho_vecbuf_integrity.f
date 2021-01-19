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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
C This file contains routines for integrity checks of the Cholesky
C vector buffer. For debugging purposes.
      Subroutine Cho_VecBuf_EnableIntegrityCheck(irc)
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Enable integrity check of buffer: allocate and store norm and sum
C     of each vector in the buffer.
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      Implicit None
      Integer irc
#include "cholesky.fh"
#include "chovecbuf.fh"
#include "WrkSpc.fh"
#include "choprint.fh"

      Real*8   dDot_, Cho_dSumElm
      external dDot_, Cho_dSumElm

      Integer ip
      Integer iSym
      Integer jVec
      Integer jRed
      Integer ipV

      ! Set return code
      irc=0

      ! Not implemented in internal run mode (since the buffered vectors
      ! may change length and hence norm and sum, causing too much book
      ! keeping activity for a debug feature).
      If (RUN_MODE.NE.RUN_EXTERNAL) Return

      ! Return if no buffer is allocated
      If (l_ChVBuf.lt.1) Return

      ! Return if already enabled
      If (l_ChVBfI.gt.0) Return

      ! Check that nDimRS is allocated
      If (.NOT.Allocated(nDimRS)) Then
         irc=1
         Return
      End If

      ! Allocate and store norm and sum of each vector in the buffer
      l_ChVBfI=0
      Do iSym=1,nSym
         l_ChVBfI_Sym(iSym)=2*nVec_in_Buf(iSym)
         l_ChVBfI=l_ChVBfI+l_ChVBfI_Sym(iSym)
      End Do
      If (l_ChVBfI.gt.0) Then
         Call Cho_Mem('ChVBfI','Allo','Real',ip_ChVBfI,l_ChVBfI)
         ip=ip_ChVBfI
         Do iSym=1,nSym
            ip_ChVBfI_Sym(iSym)=ip
            ip=ip+l_ChVBfI_Sym(iSym)
         End Do
         Do iSym=1,nSym
            ipV=ip_ChVBuf_Sym(iSym)
            ip=ip_ChvBfI_Sym(iSym)
            Do jVec=1,nVec_in_Buf(iSym)
               jRed=InfVec(jVec,2,iSym)
               Work(ip)=sqrt(dDot_(nDimRS(iSym,jRed),
     &                            Work(ipV),1,Work(ipV),1))
               Work(ip+1)=Cho_dSumElm(Work(ipV),nDimRS(iSym,jRed))
               ipV=ipV+nDimRS(iSym,jRed)
               ip=ip+2
            End Do
         End Do
         If (iPrint.gt.2) Then
            Call Cho_VecBuf_PrtRef('@NABLE')
         End If
         Write(LuPri,'(A)')
     &   'Cholesky vector buffer integrity checks enabled'
      Else
         l_ChVBfI=0
         ip_ChVBfI=0
         Call Cho_iZero(l_ChVBfI_Sym,nSym)
         Call Cho_iZero(ip_ChVBfI_Sym,nSym)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine Cho_VecBuf_Check()
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Check buffer integrity and stop if corrupted.
C
      Implicit None
#include "cholesky.fh"
      Real*8 Tol
      Logical Verbose
      Character*1 Txt
      Integer irc

      Tol=1.0d-12
      Verbose=.False.
      Txt=' '
      Call Cho_VecBuf_CheckIntegrity(Tol,Verbose,Txt,irc)
      If (irc.ne.0) Then
         Write(LuPri,'(A,I3)')
     &   'Cho_VecBuf_Check: buffer integrity check returned code',irc
         Call Cho_Quit('Cholesky vector buffer corrupted',104)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine Cho_VecBuf_CheckIntegrity(Tol,Verbose,Txt,irc)
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Check integrity of Cholesky vector buffer and print result.
C     Tol is the tolerance to use for comparing norm and sum.
C     Txt is printed along with the check message if Verbose=.True.
C     (if Verbose=.False. nothing is printed).
C     Return code:
C        irc=0: buffer OK
C        irc=1: buffer corrupted
C
C     A simpler interface is given by Subroutine Cho_VecBuf_Check.
C
      Implicit None
      Real*8  Tol
      Logical Verbose
      Character*(*) Txt
      Integer irc
#include "cholesky.fh"
      Logical Cho_VecBuf_Integrity_OK

      If (Cho_VecBuf_Integrity_OK(Tol,Verbose)) Then
         If (Verbose) Then
            Write(LuPri,'(A,A)')
     &      Txt,' Cholesky vector buffer integrity checked: OK'
            Call Cho_Flush(LuPri)
         End If
         irc=0
      Else
         If (Verbose) Then
            Write(LuPri,'(A,A)')
     &      Txt,' Cholesky vector buffer integrity checked: CORRUPTED'
            Call Cho_Quit('Buffer corrupted',104)
         End If
         irc=1
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Logical Function Cho_VecBuf_Integrity_OK(Tol,Report)
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Check Cholesky vector buffer integrity: compute norm and sum of
C     vectors in the buffer and compare these values to the table
C     generated at buffer initialization.
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      Implicit None
      Real*8  Tol
      Logical Report
#include "cholesky.fh"
#include "chovecbuf.fh"
#include "WrkSpc.fh"

      Logical OK

      Real*8   dDot_, Cho_dSumElm
      external ddot_, Cho_dSumElm


      Real*8 Nrm, Sm

      Integer nErr, iSym, jVec, jRed, n, ipV

      Integer i, j
      Real*8  RefNrm
      Real*8  RefSm
      RefNrm(i,j)=Work(ip_ChVBfI_Sym(j)+2*(i-1))
      RefSm(i,j)=Work(ip_ChVBfI_Sym(j)+2*(i-1)+1)

      nErr=0
      If (l_ChVBuf.gt.0 .and. l_ChVBfI.gt.0 .and.
     &    Allocated(nDimRS)) Then
         Do iSym=1,nSym
            If (nVec_in_Buf(iSym).gt.0 .and.
     &          l_ChVBfI_Sym(iSym).gt.0) Then
               ipV=ip_ChVBuf_Sym(iSym)
               Do jVec=1,nVec_in_Buf(iSym)
                  jRed=InfVec(jVec,2,iSym)
                  n=nDimRS(iSym,jRed)
                  Nrm=sqrt(dDot_(n,Work(ipV),1,Work(ipV),1))
                  Sm=Cho_dSumElm(Work(ipV),n)
                  OK=abs(Nrm-RefNrm(jVec,iSym)).lt.Tol .and.
     &               abs(Sm-RefSm(jVec,iSym)).lt.Tol
                  If (.not.OK) Then
                     nErr=nErr+1
                     If (Report) Then
                        Write(LuPri,'(A,I7,A,I2,A,I9)')
     &                  'Buffer corrupted: vector',jVec,' sym.',iSym,
     &                  ' dim.',n
                        Write(LuPri,'(3X,1P,3(A,D25.16))')
     &                  'Norm=',Nrm,' Reference=',RefNrm(jVec,iSym),
     &                  ' Diff=',Nrm-RefNrm(jVec,iSym)
                        Write(LuPri,'(3X,1P,3(A,D25.16))')
     &                  'Sum= ',Sm,' Reference=',RefSm(jVec,iSym),
     &                  ' Diff=',Sm-RefSm(jVec,iSym)
                     End If
                  End If
                  ipV=ipV+n
               End Do
            End If
         End Do
      End If
      If (Report) Then
         If (nErr.ne.0) Then
            Write(LuPri,'(A,I7,A,1P,D25.16)')
     &      'Buffer corrupted for ',nErr,' vectors. Tolerance=',Tol
         Else
            Write(LuPri,'(A,1P,D25.16)')
     &      'Buffer integrity OK. Tolerance=',Tol
         End If
      End If
      Cho_VecBuf_Integrity_OK=nErr.eq.0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine Cho_VecBuf_CompareNormAndSum(n,nVec,Vec,J1,iSym,irc)
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Compare norm and sum of vectors against values stored in buffer
C     for Cholesky vectors J1,J1+1,J1+nVec-1 of symmetry iSym stored in
C     array Vec(n,nVec). Comparison is only made for those vectors that
C     are actually stored in the buffer, of course.
C
C     Return codes:
C        irc=0: no differences detected.
C        irc>0: number of vectors for which differences are detected.
C
      Implicit None
      Integer n
      Integer nVec
      Real*8  Vec(n,nVec)
      Integer J1
      Integer iSym
      Integer irc
#include "WrkSpc.fh"
#include "chovecbuf.fh"

      Real*8   dDot_, Cho_dSumElm
      external ddot_, Cho_dSumElm

      Real*8 Tol
      Parameter (Tol=1.0d-12)

      Integer J0
      Integer mVec
      Integer J

      Real*8 Nrm
      Real*8 Sm

      Integer k, l
      Real*8  RefNorm
      Real*8  RefSum
      RefNorm(k,l)=Work(ip_ChVBfI_Sym(l)+2*(k-1))
      RefSum(k,l)=Work(ip_ChVBfI_Sym(l)+2*(k-1)+1)

      irc=0
      If (l_ChvBfI.gt.0) Then
         J0=J1-1
         mVec=min(J0+nVec,nVec_in_Buf(iSym))-J0
         Do J=1,mVec
            Nrm=sqrt(dDot_(n,Vec(1,J),1,Vec(1,J),1))
            Sm=Cho_dSumElm(Vec(1,J),n)
            If (abs(RefNorm(J0+J,iSym)-Nrm).gt.Tol .or.
     &          abs(RefSum(J0+J,iSym)-Sm).gt.Tol) Then
               irc=irc+1
            End If
         End Do
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine Cho_VecBuf_PrtRef(Txt)
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Print reference norm and sum of vectors in buffer.
C     Txt is printed along with the reference values (for
C     identification).
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      Implicit None
      Character*(*) Txt
#include "WrkSpc.fh"
#include "cholesky.fh"
#include "chovecbuf.fh"


      Integer iSym, jVec, jRed, nDim

      Integer i, j
      Real*8  RefNrm
      Real*8  RefSm
      RefNrm(i,j)=Work(ip_ChVBfI_Sym(j)+2*(i-1))
      RefSm(i,j)=Work(ip_ChVBfI_Sym(j)+2*(i-1)+1)

      If (.NOT.Allocated(nDimRS)) Then
         Call Cho_Quit(
     &        'Cho_VecBuf_PrtRef: unable to print reference values',104)
      End If
      If (l_ChVBfI.gt.0) Then
         Do iSym=1,nSym
            Do jVec=1,nVec_in_Buf(iSym)
               jRed=InfVec(jVec,2,iSym)
               nDim=nDimRS(iSym,jRed)
               Write(LuPri,'(A,A,I6,A,I2,A,I9,1P,2(A,D25.16))')
     &         Txt,' Cholesky vector',jVec,
     &         ' sym.',iSym,' dim.',nDim,
     &         '  Norm=',RefNrm(jVec,iSym),' Sum=',RefSm(jVec,iSym)
            End Do
         End Do
      Else
         Write(LuPri,'(A,A)')
     &   Txt,' Cho_VecBuf_PrtRef: no reference values available!'
      End If

      End
