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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_X_RdRst(ifail)
C
C     T.B. Pedersen, july 2004.
C
C     Purpose: read decomposition info and store in common
C              block. If ifail != 0 on exit, some error occurred and,
C              most likely, some of the restart info is not
C              defined/initialized.
C
#include "implicit.fh"
#include "choorb.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*11 SecNam
      Parameter (SecNam = 'Cho_X_RdRst')

      Parameter (lScr = 8)
      Real*8  dScr(lScr)
      Integer jScr(lScr)

C     Set return code.
C     ----------------

      ifail = 0

C     Read molecular info.
C     --------------------

      iAdr = 0

      iOpt = 2
      nRd  = 4
      Call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
      nShell = jScr(2)
      nnShl  = jScr(3)
      If (jScr(2) .lt. 1) Then
         Write(6,'(A,A,I10)')
     &   SecNam,': #shells from restart file:',jScr(2)
         ifail = 1
         Go To 100
      End If
      nSP_UpLim = nShell*(nShell+1)/2
      If (jScr(3).lt.1 .or. jScr(3).gt.nSP_UpLim) Then
         Write(6,'(A,A,I10)')
     &   SecNam,': #shell pairs from restart file:',jScr(3)
         ifail = 1
         Go To 100
      End If
      If (jScr(1) .ne. nSym) Then
         Write(6,'(A,A,I10)')
     &   SecNam,': #irreps from restart file:',jScr(1)
         ifail = 1
         Go To 100
      Else
         iOpt = 2
         Call iDAFile(LuRst,iOpt,jScr,nSym,iAdr)
         Do iSym = 1,nSym
            If (jScr(iSym) .ne. nBas(iSym)) Then
               Write(6,'(A,A,I2,A,I10)')
     &         SecNam,': #basis functions in sym.',iSym,
     &         ' from restart file:',jScr(iSym)
               ifail = 2
               Go To 100
            End If
         End Do
      End If

C     Read decomposition configuration info.
C     --------------------------------------

      iOpt = 2
      nRd  = 2
      Call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
      If (jScr(1) .EQ. 0) Then
         XScDiag = .false.
      Else If (jScr(1) .EQ. 1) Then
         XScDiag = .true.
      Else
         WRITE(6,'(A,A,I10)')
     &   SECNAM,': integer flag for screening not recognized:',jScr(1)
         ifail = 2
         Go To 100
      End If
      If (jScr(2).gt.0 .and. jScr(2).lt.3) Then
         XCho_AdrVec = jScr(2)
      Else
         WRITE(6,'(A,A,I10)')
     &   SECNAM,': vector file address mode not recognized:',jScr(2)
         ifail = 3
         Go To 100
      End If
      If (XCho_AdrVec .ne. Cho_AdrVec) Then
         WRITE(6,'(A,A,I10)')
     &   SECNAM,': vector file address mode from restart file:',
     &   XCho_AdrVec
         WRITE(6,'(A,A,I10)')
     &   SECNAM,': vector file address mode from runfile     :',
     &   Cho_AdrVec
         ifail = 3
         Go To 100
      End If

      iOpt = 2
      nRd  = 8
      Call dDAFile(LuRst,iOpt,dScr,nRd,iAdr)
      XThrCom  = dScr(1)
      XThrDiag = dScr(2)
      XDamp(1) = dScr(3)
      XDamp(2) = dScr(4)
      XSpan    = dScr(5)
      XThrNeg  = dScr(6)
      XWarNeg  = dScr(7)
      XTooNeg  = dScr(8)
      ThrCom   = XThrCom
      ThrDiag  = XThrDiag
      Damp(1)  = XDamp(1)
      Damp(2)  = XDamp(2)
      Span     = XSpan
      ThrNeg   = XThrNeg
      WarNeg   = XWarNeg
      TooNeg   = XTooNeg

C     Allocate InfVec array.
C     ----------------------

      l_InfVec = MaxVec*InfVec_N2*nSym
      Call GetMem('InfVec','Allo','Inte',ip_InfVec,l_InfVec)

C     Allocate and initialize (read) InfRed array.
C     --------------------------------------------

      iOpt = 2
      nRd  = 1
      Call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
      MaxRed = jScr(1)
      XnPass = MaxRed
      IF (MaxRed .lt. 1) Then
         Write(6,'(A,A,I10)')
     &   SecNam,': #reduced sets from restart file:',MaxRed
         ifail = 4
         Go To 100
      Else
         l_InfRed = MaxRed
         Call GetMem('InfRed','Allo','Inte',ip_InfRed,l_InfRed)
         iOpt = 2
         Call iDAFile(LuRst,iOpt,iWork(ip_InfRed),l_InfRed,iAdr)
         If (iWork(ip_InfRed) .ne. 0) Then
            Write(6,'(A,A,I10)')
     &      SecNam,': disk address of 1st reduced set:',iWork(ip_InfRed)
            ifail = 5
            Go To 100
         End If
      End If

C     Read InfVec array.
C     ------------------

      Do iSym = 1,nSym
         iOpt = 2
         nRd  = 1
         Call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
         If (jScr(1) .ne. NumCho(iSym)) Then
            Write(6,'(A,A,I2,A,I10)')
     &      SecNam,': #Cholesky vectors (sym.',iSym,'): ',NumCho(iSym)
            Write(6,'(A,A,I10)')
     &      SecNam,': ....and from restart file: ',jScr(iSym)
            ifail = 6
            Go To 100
         Else
            If (NumCho(iSym) .lt. 1) Then
               kOff = ip_InfVec + MaxVec*InfVec_N2*(iSym-1)
               Call Cho_iZero(iWork(kOff),MaxVec*InfVec_N2)
            Else
               Do j = 1,InfVec_N2
                  iOpt = 2
                  kOff = ip_InfVec + MaxVec*InfVec_N2*(iSym-1)
     &                 + MaxVec*(j-1)
                  Call iDAFile(LuRst,iOpt,iWork(kOff),NumCho(iSym),iAdr)
                  nRest = MaxVec - NumCho(iSym)
                  If (nRest .gt. 0) Then
                     kOff = ip_InfVec + MaxVec*InfVec_N2*(iSym-1)
     &                    + MaxVec*(j-1) + NumCho(iSym)
                     Call Cho_iZero(iWork(kOff),nRest)
                  End If
               End Do
            End If
         End If
      End Do

C     Return.
C     -------

  100 If (ifail .ne. 0) Then  ! failures jump to this point
         Write(6,'(A,A)')
     &   SecNam,': refusing to read more restart info!'
      End If

      End
