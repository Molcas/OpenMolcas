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
      SubRoutine CD_Diag(CD_Vec,
     &                   Restart,Converged,Thr,
     &                   ThrNeg,ThrFail,
     &                   DiaInp,Diag,Buf,
     &                   nDim,lBuf,
     &                   ErrStat,NumCho,
     &                   irc)
C
C     Thomas Bondo Pedersen, October 2004.
C
C     Purpose: set up diagonal for general Cholesky decomposition.
C
C     NB!! ThrNeg and ThrFail are supposed to be negative with
C          ThrFail < ThrNeg. Diagonals less than ThrNeg are zeroed,
C          while diagonals less than ThrFail is taken as a sign of
C          a non-positive definite matrix (i.e., decomposition will
C          fail).
C
C     Error codes, irc:
C        0 : all OK
C      201 : inconsistent input: Restart but NumCho < 0.
C            (NumCho >= 0 is allowed with Restart.)
C      202 : insufficient buffer size, lBuf.
C      203 : too negative diagonal element found (i.e., matrix
C            is non-positive definit).
C
      Implicit Real*8 (a-h,o-z)
      External  CD_Vec ! external routine for vectors
      Logical   Restart, Converged
      Dimension DiaInp(nDim), Diag(nDim), Buf(lBuf)
      Dimension ErrStat(3)

      Character*7 SecNam
      Parameter (SecNam = 'CD_Diag')

      Call qEnter(SecNam)

C     Set variables.
C     --------------

      irc = 0
      If (nDim .lt. 1) Then
         Converged = .true. ! in a sense, at least
         Go To 1 ! exit (nothing to do)
      Else
         Converged = .false.
      End If

C     Set up diagonal.
C     ----------------

      Call dCopy_(nDim,DiaInp,1,Diag,1)

      If (Restart) Then ! subtract previous vectors
         If (NumCho .gt. 0) Then

            nVec = min(NumCho,lBuf/nDim)
            If (nVec .lt. 1) Then
               irc = 202
               Go To 1 ! exit (insufficient buffer size)
            Else
               nBatch = (NumCho - 1)/nVec + 1
            End If

            Do iBatch = 1,nBatch

               If (iBatch .eq. nBatch) Then
                  NumV = NumCho - nVec*(nBatch - 1)
               Else
                  NumV = nVec
               End If

               iVec1 = nVec*(iBatch - 1) + 1
               iOpt  = 2
               Call CD_Vec(iVec1,NumV,Buf,lBuf,nDim,iOpt)

               Do jVec = 1,NumV
                  kOff = nDim*(jVec - 1)
                  Do i = 1,nDim
                     ij = kOff + i
                     Diag(i) = Diag(i) - Buf(ij)*Buf(ij)
                  End Do
               End Do

            End Do

         Else If (NumCho .lt. 0) Then

            irc = 201
            Go To 1 ! exit (inconsistent input)

         End If
      End If

C     Zero negative diagonals (fail if too negative),
C     get error statistics (min, max, and rms error),
C     check convergence based on max diagonal.
C     -----------------------------------------------

      If (Diag(1) .lt. ThrNeg) Then
         If (Diag(1) .lt. ThrFail) Then
            irc = 203
            Go To 1 ! exit (too negative diagonal)
         Else
            Diag(1) = 0.0d0
         End If
      End If
      ErrStat(1) = Diag(1)
      ErrStat(2) = Diag(1)
      ErrStat(3) = Diag(1)*Diag(1)
      Do i = 2,nDim
         If (Diag(1) .lt. ThrNeg) Then
            If (Diag(1) .lt. ThrFail) Then
               irc = 203
               Go To 1 ! exit (too negative diagonal)
            Else
               Diag(1) = 0.0d0
            End If
         End If
         ErrStat(1) = min(ErrStat(1),Diag(i))
         ErrStat(2) = max(ErrStat(2),Diag(i))
         ErrStat(3) = ErrStat(3) + Diag(i)*Diag(i)
      End Do
      xDim = dble(nDim)
      ErrStat(3) = sqrt(ErrStat(3))/xDim

      Converged = ErrStat(2) .le. Thr

    1   Continue
       Call qExit(SecNam)
      End
