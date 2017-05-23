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
#if defined (_CHO_DEBUG_)
#define _DEBUG_
#endif
      SubRoutine Cho_VecBuf_Init(Frac,lVec)
C
C     Purpose: allocate and initialize vector buffer.
C              RUN_MODE=RUN_INTERNAL: buffer used during decomposition.
C              RUN_MODE=RUN_EXTERNAL: buffer used after decomposition,
C                                     i.e. vectors are available on
C                                     disk.
C              (RUN_MODE stored in cholesky.fh)
C
      Implicit None
      Real*8  Frac
      Integer lVec(*)
#include "cholesky.fh"
#include "chovecbuf.fh"

      Character*15 SecNam
      Parameter (SecNam = 'Cho_VecBuf_Init')

      Logical LocDbg
#if defined (_DEBUG_)
      Parameter (LocDbg = .true.)
#else
      Parameter (LocDbg = .false.)
#endif

      Character*2 Unt
      Integer ip_Max, l_Max, MF
      Real*8  xMF

      If (LocDbg) Then
         Call Cho_Mem('GetMax','GetM','Real',ip_Max,l_Max)
         Write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
         Write(Lupri,*) 'Memory fraction requested for buffer: ',Frac
         Call Cho_Word2Byte(l_Max,8,xMF,Unt)
         Write(Lupri,*) 'Memory available: ',l_Max,' = ',xMF,Unt
         MF = INT(Frac*DBLE(l_Max))
         Call Cho_Word2Byte(MF,8,xMF,Unt)
         Write(Lupri,*) 'Memory fraction : ',MF,' = ',xMF,Unt
         Call Cho_Flush(Lupri)
      End If

      l_ChVBfI=0
      ip_ChVBfI=0
      Call iZero(l_ChVBfI_Sym,nSym)
      Call iZero(ip_ChVBfI_Sym,nSym)

      If (RUN_MODE .eq. RUN_INTERNAL) Then
         Call Cho_VecBuf_Init_I(Frac,lVec,LocDbg)
      Else If (RUN_MODE .eq. RUN_EXTERNAL) Then
         Call Cho_VecBuf_Init_X(Frac,LocDbg)
      Else
         Call Cho_Quit('RUN_MODE error in '//SecNam,103)
      End If

      If (LocDbg) Then
         Write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
         Call Cho_Flush(Lupri)
      End If

      End
      SubRoutine Cho_VecBuf_Init_I(Frac,lVec,LocDbg)
C
C     Purpose: allocate and initialize vector buffer.
C              (Internal run mode.)
C
      Implicit None
      Real*8  Frac
      Integer lVec(*)
      Logical LocDbg
#include "cholesky.fh"
#include "chovecbuf.fh"

      Character*17 SecNam
      Parameter (SecNam = 'Cho_VecBuf_Init_I')

      Character*2 Unt
      Real*8 xMemMax(8), x

      Integer lVecTot
      Integer ip_Max, l_Max
      Integer iSym, MemEach, MemLeft

      Logical Enough

      Integer  Cho_iSumElm
      External Cho_iSumElm

      If (LocDbg) Then
         Write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
         Write(Lupri,*) 'Memory fraction requested for buffer: ',Frac
         Write(Lupri,'(A,I8)')  'nSym: ',nSym
         Write(Lupri,'(A,8I8)') 'lVec: ',(lVec(iSym),iSym=1,nSym)
         Call Cho_Flush(Lupri)
      End If

      If (nSym.lt.1 .or. nSym.gt.8) Then
         Call Cho_Quit('nSym out of bounds in '//SecNam,102)
      End If

      x = DBLE(MaxVec)
      lVecTot = lVec(1)
      xMemMax(1) = DBLE(lVec(1))*x
      Do iSym = 2,nSym
         lVecTot = max(lVecTot,lVec(iSym))
         xMemMax(iSym) = DBLE(lVec(iSym))*x
      End Do

      If (Frac.le.0.0d0 .or. Frac.gt.1.0d0 .or. lVecTot.lt.1) Then
         l_ChVBuf  = 0
         ip_ChVBuf = 0
         Call Cho_iZero(ip_ChVBuf_Sym,nSym)
         Call Cho_iZero(l_ChVBuf_Sym,nSym)
      Else
         Call Cho_Mem('GetMax','GetM','Real',ip_Max,l_Max)
         l_ChVBuf = INT(Frac*DBLE(l_Max))
         If (l_ChVBuf.lt.nSym .or. l_ChVBuf.lt.lVecTot) Then
            l_ChVBuf  = 0
            ip_ChVBuf = 0
            Call Cho_iZero(ip_ChVBuf_Sym,nSym)
            Call Cho_iZero(l_ChVBuf_Sym,nSym)
         Else
            MemEach = l_ChVBuf/nSym
            Enough = MemEach .gt. lVec(1)
            Do iSym = 2,nSym
               Enough = Enough .and. MemEach.gt.lVec(iSym)
            End Do
            If (.not. Enough) Then ! whole buffer for sym. 1
               l_ChVBuf_Sym(1) = l_ChVBuf
               Do iSym = 2,nSym
                  l_ChVBuf_Sym(iSym) = 0
               End Do
            Else
               MemLeft = l_ChVBuf - nSym*MemEach
               l_ChVBuf_Sym(1) = MemEach + MemLeft
               If (DBLE(l_ChVBuf_Sym(1)) .gt. xMemMax(1)) Then
                  l_ChVBuf_Sym(1) = INT(xMemMax(1))
               End If
               Do iSym = 2,nSym
                  l_ChVBuf_Sym(iSym) = MemEach
                  If (DBLE(l_ChVBuf_Sym(iSym)) .gt. xMemMax(iSym)) Then
                     l_ChVBuf_Sym(iSym) = INT(xMemMax(iSym))
                  End If
               End Do
            End If
            l_ChVBuf = Cho_iSumElm(l_ChVBuf_Sym,nSym)
            Call Cho_Mem('CHVBUF','Allo','Real',ip_ChVBuf,l_ChVBuf)
            ip_ChVBuf_Sym(1) = ip_ChVBuf
            Do iSym = 2,nSym
               ip_ChVBuf_Sym(iSym) = ip_ChVBuf_Sym(iSym-1)
     &                             + l_ChVBuf_Sym(iSym-1)
            End Do
         End If
      End If

      Call Cho_iZero(nVec_in_Buf,nSym)

      If (LocDbg) Then
         Call Cho_Word2Byte(l_ChVBuf,8,x,Unt)
         Write(Lupri,*) 'Memory allocated for buffer: ',l_ChVBuf,
     &                  '(',x,Unt,') at ',ip_ChVBuf
         Write(Lupri,'(A,8I8)') 'l_ChVBuf_Sym : ',
     &                          (l_ChVBuf_Sym(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I8)') 'ip_ChVBuf_Sym: ',
     &                          (ip_ChVBuf_Sym(iSym),iSym=1,nSym)
         Write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
         Call Cho_Flush(Lupri)
      End If

      End
      SubRoutine Cho_VecBuf_Init_X(Frac,LocDbg)
C
C     Purpose: allocate and initialize vector buffer.
C              (External run mode.)
C
      Implicit None
      Real*8  Frac
      Logical LocDbg
#include "cholesky.fh"
#include "chovecbuf.fh"

      Character*17 SecNam
      Parameter (SecNam = 'Cho_VecBuf_Init_X')

      Integer  Cho_iSumElm
      External Cho_iSumElm

      Logical DoRead
      Integer i, iSym, ip_Max, l_Max, Left, jNum, iRedC, mUsed

      Integer lScr
      Parameter (lScr = 1)
      Real*8 Scr(lScr)

      Integer nErr
      Real*8 Scr_Check, Tol, Diff
      Parameter (Scr_Check = 1.23456789d0, Tol = 1.0d-15)

      Character*2 Unt
      Real*8 Byte

      If (LocDbg) Then
         Do i = 1,lScr
            Scr(i) = Scr_Check
         End Do
         Write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
         Write(Lupri,*) 'Memory fraction requested for buffer: ',Frac
         Write(Lupri,'(A,I8)')  'nSym: ',nSym
         Call Cho_Flush(Lupri)
      End If

      If (nSym.lt.1 .or. nSym.gt.8) Then
         Call Cho_Quit('nSym out of bounds in '//SecNam,102)
      End If

      If (Frac.le.0.0d0 .or. Frac.gt.1.0d0) Then
         l_ChVBuf  = 0
         ip_ChVBuf = 0
         Call Cho_iZero(l_ChvBuf_Sym,nSym)
         Call Cho_iZero(ip_ChvBuf_Sym,nSym)
      Else
         Call Cho_Mem('GetMax','GetM','Real',ip_Max,l_Max)
         Left = INT(Frac*DBLE(l_Max))
         iRedC = -1
         DoRead = .false.
         Do iSym = 1,nSym
            jNum = 0
            mUsed = 0
            Call Cho_VecRd1(Scr,Left,1,NumCho(iSym),iSym,
     &                      jNum,iRedC,mUsed,DoRead)
            Left = Left - mUsed
            l_ChVBuf_Sym(iSym) = mUsed
         End Do
         l_ChVBuf = Cho_iSumElm(l_ChVBuf_Sym,nSym)
         If (l_ChVBuf .lt. 1) Then
            l_ChVBuf  = 0
            ip_ChVBuf = 0
            Call Cho_iZero(l_ChvBuf_Sym,nSym)
            Call Cho_iZero(ip_ChvBuf_Sym,nSym)
         Else
            Call Cho_Mem('CHVBUF','Allo','Real',ip_ChVBuf,l_ChVBuf)
            ip_ChVBuf_Sym(1) = ip_ChVBuf
            Do iSym = 2,nSym
               ip_ChVBuf_Sym(iSym) = ip_ChVBuf_Sym(iSym-1)
     &                             + l_ChVBuf_Sym(iSym-1)
            End Do
         End If
      End If

      If (LocDbg) Then
         nErr = 0
         Do i = 1,lScr
            Diff = Scr(i) - Scr_Check
            If (ABS(Diff) .gt. Tol) Then
               nErr = nErr + 1
            End If
         End Do
         If (nErr .ne. 0) Then
            Call Cho_Quit('Memory boundary error in '//SecNam,101)
         End If
         Call Cho_Word2Byte(l_ChVBuf,8,Byte,Unt)
         Write(Lupri,*) 'Memory allocated for buffer: ',l_ChVBuf,
     &                  '(',Byte,Unt,')  at ',ip_ChVBuf
         Write(Lupri,'(A,8I8)') 'l_ChVBuf_Sym : ',
     &                          (l_ChVBuf_Sym(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I8)') 'ip_ChVBuf_Sym: ',
     &                          (ip_ChVBuf_Sym(iSym),iSym=1,nSym)
         Write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
         Call Cho_Flush(Lupri)
      End If

      End
