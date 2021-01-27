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
* Copyright (C) 2016, Thomas Bondo Pedersen                            *
************************************************************************
#if defined (_CHO_DEBUGPRINT_)
#define _DEBUGPRINT_
#endif
      SubRoutine Cho_VecBuf_Ini2()
C
C     Thomas Bondo Pedersen, June 2006.
C
C     Purpose: read vectors from disk into buffer.
C
      use ChoVecBuf
      Implicit None
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam = 'Cho_VecBuf_Ini2')

      Logical LocDbg
#if defined (_DEBUGPRINT_)
      Parameter (LocDbg = .true.)
#else
      Parameter (LocDbg = .false.)
#endif
      Logical DoRead

      Integer iV1, iV2, iSym, nRead, iRedC
      Integer mUsed(8)
      Integer irc

C     Check if buffer is allocated.
C     Check if there are any vectors.
C     -------------------------------

      If (.NOT.Allocated(CHVBUF)) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returning immediately: ',
     &                     'No buffer allocated!'
         End If
         Return
      End If
      If (NumChT .lt. 1) Then
         Write(Lupri,*) SecNam,': returning immediately: ',
     &                  'Buffer allocated, but no vectors!?!?'
         Return
      End If

C     Read vectors.
C     -------------

      DoRead = .true.
      iRedC = -1
      Do iSym = 1,nSym
         iV1 = 1
         iV2 = NumCho(iSym)
         nRead = 0
         mUsed(iSym) = 0
         Call Cho_VecRd1(CHVBUF(ip_ChVBuf_Sym(iSym)),
     &                   l_ChVBuf_Sym(iSym),
     &                   iV1,iV2,iSym,nRead,iRedC,mUsed(iSym),DoRead)
         nVec_in_Buf(iSym) = nRead
      End Do

C     Debug:
C     Enable integrity checks.
C     Print info.
C     ------------------------

      If (LocDbg) Then
         Call Cho_VecBuf_EnableIntegrityCheck(irc)
         If (irc.ne.0) Then
            Write(LuPri,'(A,I9)')
     &      SecNam,': Cho_VecBuf_EnableIntegrityCheck returned code',irc
            Call Cho_Quit(SecNam//': integrity check init failed',104)
         Else
            Write(LuPri,'(A,A)')
     &      SecNam,': buffer integrity check enabled'
         End If
         Write(Lupri,'(A,A,8I10)') SecNam,'(exit): NumCho:',
     &                             (NumCho(iSym),iSym=1,nSym)
         Write(Lupri,'(A,A,8I10)') SecNam,'(exit): nVec_in_Buf:',
     &                             (nVec_in_Buf(iSym),iSym=1,nSym)
         Write(Lupri,'(A,A,8I10)') SecNam,'(exit): buffer allocated:',
     &                             (l_ChVBuf_Sym(iSym),iSym=1,nSym)
         Write(Lupri,'(A,A,8I10)') SecNam,'(exit): memory used:',
     &                             (mUsed(iSym),iSym=1,nSym)
         Call Cho_Flush(Lupri)
      End If

      End
