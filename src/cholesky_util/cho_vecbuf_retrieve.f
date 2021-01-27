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
* Copyright (C) 2006, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_VecBuf_Retrieve(Vec,lVec,jVec1,iVec2,iSym,
     &                               jNum,iRedC,mUsed)
C
C     Thomas Bondo Pedersen, June 2006.
C
C     Purpose: copy as many vectors as possible from buffer to array
C              Vec, starting at vector jVec1 and copying at most until
C              vector iVec2. On exit, jNum is the number of vectors
C              copied. On entry as well as on exit, iRedC identifies the
C              reduced set stored at location 3 (use "-1" if none or
C              unknown). On exit, mUsed is the actual amount of memory
C              used (in array Vec).
C
C     NOTE: if no vectors can be copied, jNum=0 and mUsed=0 are returned
C           but execution is NOT stopped here!!!
C
C     NOTE: it is assumed that the vectors are stored in their
C           respective reduced sets (thus, should only be used with
C           RUN_MODE = RUN_EXTERNAL).
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      use ChoVecBuf
      Implicit Real*8 (a-h,o-z)
      Real*8 Vec(lVec)
#include "cholesky.fh"
#include "WrkSpc.fh"

      external ddot_

      Character*19 SecNam
      Parameter (SecNam = 'Cho_VecBuf_Retrieve')

      Logical Full
      Logical LocDbg
#ifdef _DEBUGPRINT_
      Parameter (LocDbg = .true.)
#else
      Parameter (LocDbg = .false.)
#endif

C     Initialize.
C     -----------

      jNum = 0
      mUsed = 0

C     Check that a buffer has been allocated and contains vectors in the
C     requested range.
C     ------------------------------------------------------------------

      If (l_ChvBuf_Sym(iSym) .lt. 1) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returning immediately. ',
     &                     'No buffer allocated.'
         End If
         Return
      End If
      If (l_ChvBfI_Sym(iSym).gt.0) Then
         Call Cho_VecBuf_Check()
      End If
      If (jVec1 .gt. nVec_in_Buf(iSym)) Then
         If (LocDbg) Then
            Write(Lupri,*) SecNam,': returning immediately. ',
     &                     'jVec1 = ',jVec1,'  >  nVec_in_Buf = ',
     &                     nVec_in_Buf(iSym),' (sym. ',iSym,')'
         End If
         Return
      End If

C     Count how many vectors can be copied.
C     -------------------------------------

      lTot = 0
      Full = lTot .ge. lVec
      jVec = jVec1 - 1
      iV2  = min(nVec_in_Buf(iSym),iVec2)
      If (.NOT.Allocated(nDimRS)) Then
         iLoc = 3
         Do While (jVec.lt.iV2 .and. .not.Full)
            jVec = jVec + 1
            jRed = InfVec(jVec,2,iSym)
            If (jRed .ne. iRedC) Then
               irc = 0
               Call Cho_X_SetRed(irc,iLoc,jRed)
               If (irc .ne. 0) Then
                  Write(Lupri,*) SecNam,': Cho_X_SetRed returned ',irc
                  Call Cho_Quit('Error in '//SecNam,104)
               End If
               iRedC = jRed
            End If
            lTot = lTot + nnBstR(iSym,iLoc)
            If (lTot .gt. lVec) Then
               jVec = jVec - 1
               lTot = lTot - nnBstR(iSym,iLoc)
               Full = .true.
            Else
               jNum = jNum + 1
            End If
         End Do
      Else
         Do While (jVec.lt.iV2 .and. .not.Full)
            jVec = jVec + 1
            jRed = InfVec(jVec,2,iSym)
            lTot = lTot + nDimRS(iSym,jRed)
            If (lTot .gt. lVec) Then
               jVec = jVec - 1
               lTot = lTot - nDimRS(iSym,jRed)
               Full = .true.
            Else
               jNum = jNum + 1
            End If
         End Do
      End If

C     Copy vectors (if any).
C     ----------------------

      If (lTot .gt. 0) Then
         kB = ip_ChVBuf_Sym(iSym)
         If (jVec1 .gt. 1) Then
            If (.NOT.Allocated(nDimRS)) Then
               iLoc = 3
               Do jVec = 1,jVec1-1
                  jRed = InfVec(jVec,2,iSym)
                  If (iRedC .ne. jRed) Then
                     irc = 0
                     Call Cho_X_SetRed(irc,iLoc,jRed)
                     If (irc .ne. 0) Then
                        Write(Lupri,*) SecNam,
     &                                 ': Cho_X_SetRed returned ',irc
                        Call Cho_Quit('Error [2] in '//SecNam,104)
                     End If
                     iRedC = jRed
                  End If
                  kB = kB + nnBstR(iSym,iLoc)
               End Do
            Else
               Do jVec = 1,jVec1-1
                  jRed = InfVec(jVec,2,iSym)
                  kB = kB + nDimRS(iSym,jRed)
               End Do
            End If
         End If
         Call dCopy_(lTot,CHVBUF(kB),1,Vec,1)
         ! Check copy operation (may fail if molcas is compiled for
         ! 64 bit but linked to a 32 bit blas library)
         ! Note: check is not done unless it is enabled when the buffer
         ! is initialized.
         ! This only happens if _DEBUGPRINT_ is defined, so
         ! the following section is normally not executed.
         ! For debugging without turning on _DEBUGPRINT_ compilation:
         !   Call Cho_VecBuf_EnableIntegrityCheck(irc)
         ! somewhere in the code and for every subsequent call to this
         ! routine, the integrity is checked.
         If (l_ChVBfI_Sym(iSym).gt.0) Then
            nErr=0
            kB=1
            Do iVec=1,jNum
               jVec=jVec1+iVec-1
               jRed=InfVec(jVec,2,iSym)
               Call Cho_VecBuf_CompareNormAndSum(nDimRS(iSym,jRed),1,
     &                                           Vec(kB),jVec,iSym,irc)
               If (irc.ne.0) Then
                  nErr=nErr+1
                  Write(LuPri,'(A,I9,A,I2,A)')
     &            'Buffer copy failed for vector',jVec,' (sym.',iSym,')'
               End If
               kB=kB+nDimRS(iSym,jRed)
            End Do
            If (nErr.gt.0) Then
               Call Cho_Flush(LuPri)
               Write(LuPri,'(A,I9,A)')
     &         'Cho_VecBuf_Retrieve: buffer copy failed for',nErr,
     &         ' vectors. Going to check buffer integrity...'
               Call Cho_Flush(LuPri)
               Call Cho_VecBuf_Check()
               Write(LuPri,'(A,A)')
     &         'Buffer integrity checked: OK',
     &         ' --- error occurs in the copy operation.'
#if defined (_I8_)
               Write(LuPri,'(A,A)')
     &         'This appears to be a 64-bit version of MOLCAS.',
     &         ' Did you link to a 32-bit version of the BLAS library?'
#endif
               Call Cho_Quit(
     &         'Cho_VecBuf_Retrieve: buffer copy failed',104)
            End If
         End If
      End If

C     Set memory used.
C     ----------------

      mUsed = lTot

C     Debug: print.
C     -------------

      If (LocDbg) Then
         Write(Lupri,*)
         Write(Lupri,*) SecNam,':'
         If (jNum .lt. 1) Then
            Write(Lupri,*) 'No vectors copied!'
         Else
            Write(Lupri,*) 'Vectors ',jVec1,' to ',jVec1+jNum-1,
     &                     ' of symmetry ',iSym,' copied from buffer.'
            If (Allocated(nDimRS)) Then
               kOffV = 1
               Do iVec = 1,jNum
                  jVec = jVec1 + iVec - 1
                  jRed = InfVec(jVec,2,iSym)
                  jAdr = InfVec(jVec,3,iSym)
                  xNrm = sqrt(dDot_(nDimRS(iSym,jRed),Vec(kOffV),1,
     &                                               Vec(kOffV),1))
                  Write(Lupri,*) 'Vector:',jVec,' disk address: ',jAdr,
     &                           ' norm: ',xNrm
                  kOffV = kOffV + nDimRS(iSym,jRed)
               End Do
               nTst = kOffV - 1
               If (nTst .NE. mUsed) Then
                 Call Cho_Quit('Vector dimension error in '//SecNam,104)
               End If
            End If
         End If
         Call Cho_Flush(Lupri)
      End If

      End
