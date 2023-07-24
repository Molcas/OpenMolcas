!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************
      SubRoutine CD_InCore_1(X,n,Vec,MxVec,NumCho,Thr,ThrNeg,ThrFail,
     &                       irc)
!
!     Thomas Bondo Pedersen, October 2004.
!
!     Purpose: Cholesky decompose the n-by-n matrix X.
!              Vectors are returned in Vec array.
!
!     Return code (irc):
!
!        0 -- decomposition success.
!      101 -- negative diagonal encountered (i.e. < ThrFail)
!      102 -- number of vectors needed exceeds max. allowed (MxVec)
!
!     Note: the algorithm is designed for incomplete Cholesky
!     decomposition, i.e. for semi-definitive matrices, and thus makes
!     use of level-1 BLAS only.
!
      Implicit None
      Integer n, MxVec, NumCho, irc
      Real*8  X(n,n), Vec(n,MxVec)
      Real*8  Thr, ThrNeg, ThrFail

      Integer i, imax, j
      Integer iPass
      Real*8  Xmax, Factor, xFac, Acc

      irc = 0

      NumCho = 0
      Acc=Min(1.0D-12,Thr*1.0D-2)
      xFac = 0.0D0  ! dummy set
      Do iPass = 1,n

         If (X(1,1) .lt. ThrNeg) Then
            If (X(1,1) .lt. ThrFail) Then
               irc = 101
               Return
            Else
               Do j = 1,n
                  X(j,1) = 0.0d0
                  X(1,j) = 0.0d0
               End Do
            End If
         End If
         Xmax = 0.0D0
         imax = 0
         Do i = 1,n
            If (X(i,i) .lt. ThrNeg) Then
               If (X(i,i) .lt. ThrFail) Then
                  irc = 101
                  Return
               Else
                  Do j = 1,n
                     X(j,i) = 0.0d0
                     X(i,j) = 0.0d0
                  End Do
               End If
            End If
            If (X(i,i) .gt. Xmax+Acc) Then
               Xmax = X(i,i)
               imax = i
            End If
         End Do

         If (Xmax .le. Thr) Return ! converged
         xFac = 1.0d0/sqrt(abs(Xmax))

         If (NumCho .eq. MxVec) Then ! too many vectors
            irc = 102
            Return
         End If

         Do j = 1,NumCho
            Factor = -Vec(imax,j)
            Call dAXPY_(n,Factor,Vec(1,j),1,X(1,imax),1)
         End Do
         X(imax,imax) = Xmax

         NumCho = NumCho + 1
         Do i = 1,n
            If (X(i,i) .eq. 0.0d0) Then
               Vec(i,NumCho) = 0.0d0
            Else
               Vec(i,NumCho) = xFac*X(i,imax)
            End If
         End Do

         Do i = 1,n
            X(i,i) = X(i,i) - Vec(i,NumCho)**2
         End Do
         X(imax,imax) = 0.0d0

      End Do

      End
