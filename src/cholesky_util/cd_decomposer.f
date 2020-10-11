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
* Copyright (C) 2004,2008, Thomas Bondo Pedersen                       *
************************************************************************
      SubRoutine CD_Decomposer(CD_Col,CD_Vec,MxNumCho,
     &                         Thr,Span,MxQual,
     &                         ThrNeg,ThrFail,
     &                         Diag,Qual,Buf,
     &                         iPivot,iQual,
     &                         nDim,lBuf,
     &                         NumCho,
     &                         irc)

C
C     Thomas Bondo Pedersen, October 2004.
C     Modified to compute at most MxNumCho vectors,
C        Thomas Bondo Pedersen, January 2008.
C
C     Purpose: Cholesky decompose a matrix.
C              Stop decomposition when either
C              1) max. diag <= Thr
C              2) NumCho = MxNumCho
C
C     To use criterion 1) only (standard procedure),
C     simply set MxNumCho = nDim.
C     To use criterion 2) only,
C     simply set Thr=1.0d-20 (i.e. zero)
C
C     Note: do *not* call this routine directly;
C           use ChoDec(...) or ChoDec_MxVec instead
C           (see those routines for documentation).
C           This routine contains implicit assumptions
C           that are checked by ChoDec and ChoDec_MxVec!!!
C
C     Error codes, irc:
C        0 : all OK
C      301 : too few qualified (probably a bug)
C      302 : insufficient buffer size, lBuf
C      303 : too negative diagonal encountered
C            (matrix non-positive definite!)
C

      Implicit Real*8 (a-h,o-z)
      External CD_Col    ! external routine for matrix columns
      External CD_Vec    ! external routine for Cholesky vectors
      Dimension Diag(nDim), Qual(nDim,MxQual), Buf(lBuf)
      Integer   iPivot(nDim), iQual(MxQual)

      Character*13 SecNam
      Parameter (SecNam = 'CD_Decomposer')

      Logical Last

      irc = 0

      iPass = 0
      mPass = MxNumCho
      Do While (iPass .lt. mPass)

C        Update counter.
C        ---------------

         iPass = iPass + 1

C        Find max. diagonal.
C        -------------------

         Dmax = Diag(1)
         Do i = 2,nDim
            Dmax = max(Dmax,Diag(i))
         End Do

C        Check for convergence.
C        ----------------------

         If (Dmax.gt.Thr .and. NumCho.lt.MxNumCho) Then

C           Find largest diagonal elements > DiaMin.
C           I.e., qualify columns.
C           ========================================

            nQual  = min(MxQual,MxNumCho-NumCho)
            DiaMin = max(Dmax*Span,Thr)
            Call CD_DiaMax(Diag,nDim,iPivot,iQual,nQual,DiaMin)

            If (nQual .lt. 1) Then ! this would be a bug...
               irc = 301
               Go To 1  ! exit
            End If

C           Get qualified columns from external routine.
C           ============================================

            Call CD_Col(Qual,nDim,iQual,nQual,Buf,lBuf)

C           Subtract previous vectors (if any).
C           ===================================

            If (NumCho .gt. 0) Then

               MinBuf = nDim + nQual
               nVec   = min(NumCho,lBuf/MinBuf)
               If (nVec .lt. 1) Then  ! insufficient buffer size
                  irc = 302
                  Go To 1 ! exit
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
                  lVec  = nDim*NumV

                  kOffV = 1
                  kOffQ = kOffV + lVec

                  iOpt = 2
                  Call CD_Vec(iVec1,NumV,Buf(kOffV),lBuf-kOffV+1,nDim,
     &                        iOpt)

                  Do jVec = 1,NumV
                     Do i = 1,nQual
                        kOff1 = kOffQ + NumV*(i - 1) + jVec - 1
                        kOff2 = kOffV + nDim*(jVec - 1)
     &                        + iQual(i) - 1
                        Buf(kOff1) = Buf(kOff2)
                     End Do
                  End Do

                  Call DGEMM_('N','N',nDim,nQual,NumV,
     &                       -1.0d0,Buf(kOffV),nDim,Buf(kOffQ),NumV,
     &                       1.0d0,Qual(1,1),nDim)

               End Do

            End If

C           Decompose.
C           ==========

            MxVec  = min(nQual,lBuf/nDim)
            iDump  = 0
            iChoMx = nQual
            iCho   = 0
            Do While (iCho .lt. iChoMx)

C              Find max. among qualified.
C              --------------------------

               Dx = Diag(iQual(1))
               ix = 1
               Do i = 2,nQual
                  If (Diag(iQual(i)) .gt. Dx) Then
                     Dx = Diag(iQual(i))
                     ix = i
                  End If
               End Do

               Last = Dx.lt.DiaMin .or. Dx.le.Thr
               If (.not. Last) Then

C                 Calculate new vector.
C                 ---------------------

                  Factor = 1.0d0/sqrt(Dx)
                  Do i = 1,nDim
                     If (Diag(i) .eq. 0.0d0) Then
                        Qual(i,ix) = 0.0d0
                     Else
                        Qual(i,ix) = Factor*Qual(i,ix)
                     End If
                  End Do

C                 Update diagonal and find new max.
C                 ---------------------------------

                  Diag(1) = Diag(1) - Qual(1,ix)*Qual(1,ix)
                  xm = Diag(1)
                  Do i = 2,nDim
                     Diag(i) = Diag(i) - Qual(i,ix)*Qual(i,ix)
                     xm = max(xm,Diag(i))
                  End Do

C                 Zero treated diagonal and find new DiaMin.
C                 ------------------------------------------

                  Diag(iQual(ix)) = 0.0d0
                  DiaMin = max(xm*Span,Thr)

C                 Zero negative diagonals (quit if too negative).
C                 -----------------------------------------------

                  Do i = 1,nDim
                     If (Diag(i) .lt. ThrNeg) Then
                        If (Diag(i) .lt. ThrFail) Then
                           irc = 303
                           Go To 1 ! exit (too negative diagonal)
                        Else
                           Diag(i) = 0.0d0
                        End If
                     End If
                  End Do

C                 Subtract this vector from qualified columns.
C                 --------------------------------------------

                  Do i = 1,nQual
                     If (Diag(iQual(i)) .ne. 0.0d0) Then
                        Factor = -Qual(iQual(i),ix)
                        Call dAXPY_(nDim,Factor,
     &                             Qual(1,ix),1,Qual(1,i),1)
                     End If
                  End Do

C                 Store vector in buffer.
C                 -----------------------

                  kOff = nDim*iDump + 1
                  Call dCopy_(nDim,Qual(1,ix),1,Buf(kOff),1)

C                 Update counter.
C                 ---------------

                  iDump = iDump + 1

               End If

C              Dump vectors to external routine CD_Vec.
C              ----------------------------------------

               If (Last .or. iDump.eq.MxVec) Then
                  If (iDump .gt. 0) Then
                     iVec1 = NumCho + 1
                     iOpt = 1
                     Call CD_Vec(iVec1,iDump,Buf,lBuf,nDim,iOpt)
                     NumCho = NumCho + iDump
                  End If
                  If (Last) Then
                     iCho = iChoMx + 1 ! break iCho loop (next pass)
                  Else
                     iDump = 0
                     iCho  = iCho + 1
                  End If
               End If

            End Do

         Else ! converged; break while loop

            iPass = mPass + 1

         End If

      End Do

    1 Continue
      End
