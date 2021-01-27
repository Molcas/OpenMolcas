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
      SubRoutine Cho_VecBuf_Maintain(irc,iRed,DoTime,DoStat)
C
C     Purpose: maintain Cholesky vector buffer:
C
C              1) reorder vectors in buffer to current reduced set
C                 storage (defined by location 2),
C
C              2) if possible, read in new vectors and store in current
C                 reduced set storage (defined by location 2).
C
C              It is assumed that all vectors in the buffer are stored
C              according to the reduced set identified by iRed.
C
C     DoTime: time as vector I/O.
C     DoStat: update statistics info (#calls to system for I/O).
C
C     Return code:  irc  = 0 : success
C                   irc != 0 : failure
C
C     Index arrays from chovecbuf.f90 modified by this routine:
C
C     NVEC_IN_BUF() -- #vectors stored in buffer in each symmetry
C
      use ChoArr, only: iScr
      use ChoSwp, only: InfVec
      use ChoVecBuf
#include "implicit.fh"
      Logical DoTime, DoStat
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character*19 SecNam
      Parameter (SecNam = 'Cho_VecBuf_Maintain')

      Logical LocDbg
*#define _DEBUGPRINT_
#if defined (_DEBUGPRINT_)
      Parameter (LocDbg = .true.)
#else
      Parameter (LocDbg = .false.)
#endif

C     Debug print.
C     ------------

      If (LocDbg) Then
         Write(Lupri,*)
         Write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
         Write(Lupri,*) 'iRed = ',iRed
         Write(Lupri,*) 'l_ChVBuf  = ',l_ChVBuf,
     &                  '   ip_ChVBuf = ',ip_ChVBuf
         Write(Lupri,'(A,8I16)') 'l_ChVBuf_Sym : ',
     &                          (l_ChVBuf_Sym(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I16)') 'ip_ChVBuf_Sym: ',
     &                          (ip_ChVBuf_Sym(iSym),iSym=1,nSym)
         Write(Lupri,'(A,8I16)') 'nVec_in_Buf  : ',
     &                          (nVec_in_Buf(iSym),iSym=1,nSym)
      End If

C     Set return code.
C     ----------------

      irc = 0

C     Return if there is no buffer to maintain.
C     -----------------------------------------

      If (l_ChVBuf .lt. 1) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returning: no buffer to maintain!'
            Write(Lupri,*) SecNam,': l_ChVBuf = ',l_ChVBuf
         End If
         Return
      End If

C     If there are no vectors yet, return.
C     ------------------------------------

      If (NumChT .lt. 1) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returning: no vectors!'
            Write(Lupri,*) SecNam,': NumChT = ',NumChT
         End If
         Return
      End If

C     Check that iScr array has been allocated.
C     -----------------------------------------

      If (.NOT.Allocated(iScr)) Then
         Write(Lupri,*) SecNam,': iScr array not allocated!'
         irc = 102
         Return
      End If

C     Start timing.
C     -------------

      If (DoTime) Call Cho_Timer(C1,W1)

C     Set index arrays for reduced set iRed at location 3.
C     ----------------------------------------------------

      If (iRed .lt. 1) Then
         nErr = 0
         Do iSym = 1,nSym
            If (nVec_in_Buf(iSym) .ne. 0) Then
               Write(Lupri,*) SecNam,': sym. block ',iSym,':'
               Write(Lupri,*) '   iRed = ',iRed,
     &                        '   #vectors in buffer: ',
     &                        nVec_in_Buf(iSym),' (should be 0)'
               nErr = nErr + 1
            End If
         End Do
         If (nErr .ne. 0) Then
            irc = 103
            Return
         End If
         iRedC = 1
      Else
         iRedC = iRed
      End If

      Call Cho_X_SetRed(irc,3,iRedC)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': Cho_X_SetRed returned ',irc
         irc = 104
         Return
      End If

C     Reordering.
C     ===========

      Do iSym = 1,nSym
         If (nnBstR(iSym,2).gt.0 .and. nVec_in_Buf(iSym).gt.0) Then

C           Check reduced set dimensions.
C           -----------------------------

            If (nnBstR(iSym,2) .gt. nnBstR(iSym,3)) Then
               Write(Lupri,*) SecNam,': dimension of reduced set at 2 ',
     &                        ' is larger than that at 3'
               Write(Lupri,*) 'Symmetry: ',iSym,'  ID of 3: ',iRedC
               Write(Lupri,*) 'Dimension of reduced set 2: ',
     &                        nnBstR(iSym,2)
               Write(Lupri,*) 'Dimension of reduced set 3: ',
     &                        nnBstR(iSym,3)
               Write(Lupri,*) 'Unable to continue...'
               irc = 104
               Return
            End If

C           Define mapping from reduced set at location 2 to that at
C           location 3.
C           --------------------------------------------------------

            Call Cho_RS2RS(iScr,SIZE(iScr),2,3,iRedC,iSym)

C           Reorder vectors.
C           ----------------

            Do iVec = 1,nVec_in_Buf(iSym)
               iOff2 = ip_ChVBuf_Sym(iSym) + nnBstR(iSym,2)*(iVec-1) - 1
               iOff3 = ip_ChVBuf_Sym(iSym) + nnBstR(iSym,3)*(iVec-1) - 1
               Do iRS2 = 1,nnBstR(iSym,2)
#if defined (_DEBUGPRINT_)
                  jRS3 = iScr(iRS2)
                  If (iRS2.lt.1 .or. iRS2.gt.SIZE(iScr)) Then
                     Write(LuPri,*) 'iRS2=',iRS2
                     Write(LuPri,*) 'SIZE(iScr)=',SIZE(iScr)
                     Call Cho_Quit('RS-2-RS map error in '//SecNam,104)
                  End If
                  If (jRS3.lt.1 .or. jRS3.gt.nnBstR(iSym,3)) Then
                     Write(LuPri,*) 'jRS3=',JRS3
                     Write(LuPri,*) 'nnBstR(iSym,3)=',nnBstR(iSym,3)
                     Call Cho_Quit('RS-2-RS map error in '//SecNam,104)
                  End If
#endif
                  Work(iOff2+iRS2) = Work(iOff3+iScr(iRS2))
               End Do
            End Do

         End If
      End Do

C     Read in more vectors.
C     =====================

      nSys = 0 ! #calls to reading routine (counter)

      Call Cho_Mem('CHVB.Read','MAX ','Real',ip_VRd,l_VRd)
      Do iSym = 1,nSym
         nDisk = NumCho(iSym) - nVec_in_Buf(iSym)
#if defined (_DEBUGPRINT_)
         If (nDisk .lt. 0) Then
            Call Cho_Quit('nDisk < 0 in '//SecNam,103)
         End If
#endif
         iMapC = -1
         If (nnBstR(iSym,2).gt.0 .and. nDisk.gt.0) Then

C           Compute how many more vectors can be stored in buffer taking
C           into account the number of vectors on disk.
C           ------------------------------------------------------------

            Left = l_ChVBuf_Sym(iSym) - nnBstR(iSym,2)*nVec_in_Buf(iSym)
            If (Left .ge. 0) Then
               nVec = min(Left/nnBstR(iSym,2),nDisk)
            Else
               Call Cho_Quit('Left < 0 in '//SecNam,103)
               nVec = 0
            End If
            iVec1 = nVec_in_Buf(iSym) + 1
            iVec2 = iVec1 + nVec - 1

C           Read and reorder vectors.
C           -------------------------

            iVec = iVec1
            Do While (iVec .le. iVec2)

C              Read vectors.
C              -------------

               nVRd  = 0
               mUsed = 0
               Call Cho_VecRd(Work(ip_VRd),l_VRd,iVec,iVec2,iSym,nVRd,
     &                        iRedC,mUsed)
               If (nVRd .lt. 1) Then
                  Call Cho_Quit('Insufficient memory for read in '
     &                          //SecNam,101)
               End If
               nSys = nSys + 1

C              Reorder the vectors and store in buffer in current
C              reduced set.
C              --------------------------------------------------

               iOff2 = ip_ChVBuf_Sym(iSym)
     &               + nnBstR(iSym,2)*nVec_in_Buf(iSym) - 1
               iOff3 = ip_VRd - 1
               Do kVec = 1,nVRd

                  jVec = iVec + kVec - 1
                  jRed = InfVec(jVec,2,iSym)
                  If (jRed .ne. iRedC) Then
                     Call Cho_X_SetRed(irc,3,jRed)
                     If (irc .ne. 0) Then
                        Write(Lupri,*)
     &                  SecNam,': Cho_X_SetRed [2] returned ',irc
                        irc = 104
                        Return
                     End If
                     iRedC = jRed
                  End If

                  If (jRed .ne. iMapC) Then
                     Call Cho_RS2RS(iScr,SIZE(iScr),2,3,jRed,iSym)
                     iMapC = jRed
                  End If

                  Do iRS2 = 1,nnBstR(iSym,2)
#if defined (_DEBUGPRINT_)
                     jRS3 = iScr(iRS2)
                     If (jRS3.lt.1 .or. jRS3.gt.nnBstR(iSym,3)) Then
                        Call Cho_Quit('RS-2-RS map error [2] in '
     &                                //SecNam,104)
                     End If
#endif
                     Work(iOff2+iRS2) = Work(iOff3+iScr(iRS2))
                  End Do

                  iOff2 = iOff2 + nnBstR(iSym,2)
                  iOff3 = iOff3 + nnBstR(iSym,3)

               End Do

C              Update counters.
C              ----------------

               iVec = iVec + nVRd
               nVec_in_Buf(iSym) = nVec_in_Buf(iSym) + nVRd

            End Do

         End If
      End Do
      Call Cho_Mem('CHVB.Read','Free','Real',ip_VRd,l_VRd)

C     Update global timing.
C     ---------------------

      If (DoStat) nSys_Call = nSys_Call + nSys

C     Update global timing.
C     ---------------------

      If (DoTime) Then
         Call Cho_Timer(C2,W2)
         tDecom(1,2) = tDecom(1,2) + C2 - C1
         tDecom(2,2) = tDecom(2,2) + W2 - W1
      End If

C     Debug print.
C     ------------

      If (LocDbg) Then
         Write(Lupri,*) 'After updating: '
         Write(Lupri,'(A,8I8)') 'nVec_in_Buf  : ',
     &                          (nVec_in_Buf(iSym),iSym=1,nSym)
         Write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
      End If

      End
