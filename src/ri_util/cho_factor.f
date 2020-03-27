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
* Copyright (C) 2007, Francesco Aquilante                              *
*               2014, Thomas Bondo Pedersen                            *
************************************************************************
*  CHO_FACTOR
*
*> @brief
*>   Evaluation of the Cholesky factor (\f$ Z \f$) of a SPD matrix (\f$ A \f$)
*> @author F. Aquilante (Jan. 2007)
*> @modified_by T.B. Pedersen (2014) Changed criterion for too negative diagonal
*>
*> @details
*> Evaluation of the Cholesky factor (\f$ Z \f$) of a SPD matrix (\f$ A \f$)
*>
*> \code
*>   For k=1,dim(A)
*>     Z(k,k) = sqrt( A(k,k) - sum_j  Z(k,j)^2 )
*>     Z(i,k) = ( A(i,k) - sum_j  Z(i,j)*Z(k,j) ) / Z(k,k)
*> \endcode
*>
*> The result is such that \f$ A \f$ is Cholesky decomposed as
*>
*> \f[  A = Z Z^\text{T} \f]
*>
*> The Cholesky factor is in general *NOT UNIQUE!!*
*> Therefore, and for stability reason, pivoting of the
*> initial matrix \f$ A \f$ would be advisable.
*>
*> @side_effects
*> \p A_k in output contains the \p kCol -th Cholesky vector.
*> In case of detected linear dependence, the \p A_k array
*> is returned as zeros!
*>
*> @note
*> Rectangular storage must be used for the \f$ Z \f$-matrix!
*>
*> @param[in,out] Diag   Updated diagonal elements of \f$ A \f$ (subtraction done by this routine)
*> @param[in,out] A_k    currently treated column of \f$ A \f$. In output contains the \p kCol -th Cholesky vector
*> @param[in]     iD_A   indices of the columns of \f$ A \f$
*> @param[in]     kCol   index of the Cholesky vector
*> @param[in]     nRow   number of rows of \f$ A \f$
*> @param[in]     Zm     in-core matrix whose columns are the Cholesky vectors
*> @param[in]     nMem   max number of columns of \p Zm kept in core
*> @param[in]     lu_Z   file unit where the \f$ Z \f$-matrix is stored
*> @param[in]     Scr    scratch space used for reading out-of-core columns of \p Zm
*> @param[in]     lScr   size of the scratch space (&ge; \p nRow or ``0`` iff in-core)
*> @param[in]     thr    threshold for linear dependence
*> @param[out]    lindep integer indicating detected linear dependence (= ``1`` iff found lin dep, else = ``0``)
************************************************************************
      SUBROUTINE CHO_FACTOR(Diag,A_k,iD_A,kCol,nRow,Zm,nMem,lu_Z,Scr,
     &                      lScr,thr,lindep)

      Implicit Real*8 (a-h,o-z)
      Integer iD_A(*), kCol, nRow, nMem, lu_Z, lScr, lindep
      Real*8  Diag(*), A_k(*), Zm(nRow,*), Scr(*)
#include "warnings.fh"
      Parameter ( one = 1.0d0, zero = 0.0d0 , thr_neg=-1.0d-8)

************************************************************************

      If (thr .lt. zero) Then
         Call WarningMessage(2,'Error in Cho_Factor')
         write(6,*)'thr must be .ge. zero'
         Call Quit(_RC_CHO_LOG_)
      EndIf

      lindep = 0
      Dmax = Diag(iD_A(kCol)) ! pivoting done by the calling routine
      xfac = one/sqrt(Abs(Dmax))

      If (kCol .le. nMem) Then

         If (Dmax.ge.thr) Then
*
*  Compute elements of the k-th Cholesky vector
* -------------------------------------------------------
C     Z(i,k) = A(i,k) - sum_j  Z(k,j)*Z(i,j)
* -------------------------------------------------------
            Do j=1,kCol-1

               fac = -Zm(iD_A(kCol),j)
               Call dAXPY_(nRow,fac,Zm(1,j),1,A_k(1),1)

            End Do

C-tbp: use thr_neg as threshold for too negative diagonal
C      It should not depend on the decomposition threshold!
C        ElseIf (Dmax.gt.zero .or. -Dmax.le.1.0d1*thr) Then
         ElseIf (Dmax.gt.thr_neg) Then

            lindep = 1
            Call Fzero(A_k(1),nRow)
            Return

         Else

            Call WarningMessage(2,'Error in Cho_Factor')
            write(6,*)'CHO_FACTOR: too-negative diagonal.'
            write(6,*)'CHO_FACTOR: current largest Diag = ',Dmax
            Call Quit(_RC_CHO_RUN_)

         EndIf
*                                                                      *
************************************************************************
*                                                                      *
      Else   ! the first nMem columns of Z are in memory
*                                                                      *
************************************************************************
*                                                                      *

         If (lScr .lt. nRow) Then
            Call WarningMessage(2,'Error in Cho_Factor')
            write(6,*)'lScr must be .ge. nRow'
            Call Quit(_RC_CHO_LOG_)
         EndIf

         If (Dmax.ge.thr) Then
*
*  Compute elements of the k-th Cholesky vector (in-core contrib.)
* -----------------------------------------------------------------
C     Z(i,k) = A(i,k) - sum_j  Z(k,j)*Z(i,j)
* -----------------------------------------------------------------
            Do j=1,nMem

               fac = -Zm(iD_A(kCol),j)
               Call dAXPY_(nRow,fac,Zm(1,j),1,A_k(1),1)

            End Do
*
*  Batch for the out-of-core previous vectors
*--------------------------------------------
            kstep = lScr/nRow

            Do kdone=nMem+1,kCol-1,kStep

               lZdone = nRow*(kdone-1)
               lZrem = nRow*(kCol-kdone)
               lZread = Min(LZrem,nRow*kStep)

               Call ddafile(lu_Z,2,Scr(1),lZread,lZdone) ! read
*
*  Compute elements of the k-th Cholesky vector (out-of-core contrib.)
* --------------------------------------------------------------------
C     Z(i,k) = A(i,k) - sum_j  Z(k,j)*Z(i,j)
* --------------------------------------------------------------------
               Do j=1,lZread/nRow
                  kj = nRow*(j-1)
                  fac = -Scr(kj+iD_A(kCol))
                  Call dAXPY_(nRow,fac,Scr(1+kj),1,A_k(1),1)
               End Do

            End Do

C-tbp: use thr_neg as threshold for too negative diagonal
C      It should not depend on the decomposition threshold!
C-tbp    ElseIf (Dmax.gt.zero .or. -Dmax.le.1.0d1*thr) Then
         ElseIf (Dmax.gt.thr_neg) Then

            lindep = 1
            Call Fzero(A_k(1),nRow)
            Return

         Else

            Call WarningMessage(2,'Error in Cho_Factor')
            write(6,*)'CHO_FACTOR: too-negative diagonal.'
            write(6,*)'CHO_FACTOR: current largest Diag = ',Dmax
            Call Quit(_RC_CHO_RUN_)

         EndIf

      EndIf
*
      A_k(iD_A(kCol)) = Dmax
*
* Scaling of the vector elements :  Z(i,k) = Z(i,k)/Z(k,k)
* --------------------------------------------------------
      call dscal_(nRow,xfac,A_k(1),1)
*
*  Explicit zeroing of the previously treated elements
* ----------------------------------------------------
      Do i=1,kCol-1
         A_k(iD_A(i)) = zero
      End Do
*
*  Update diagonal elements of the A matrix
*  ----------------------------------------------------
C     A(i,i) = A(i,i) - Z(i,k)^2    ( i > k )
*  ----------------------------------------------------
      Do i=1,nRow
         Diag(i) = Diag(i) - A_k(i)**2
      End Do
      Diag(iD_A(kCol)) = zero ! explicit zeroing of the treated diagonal

C-tbp: zero negative diagonal elements
C      Stop if too negative!
      Do i = 1,nRow
         If (Diag(i).lt.zero) Then
            If (Diag(i).le.thr_neg) Then
               Call WarningMessage(2,'Error in Cho_Factor')
               write(6,*)'CHO_FACTOR: too negative diagonal.'
               write(6,*)'CHO_FACTOR: i,Diag(i)= ',i,Diag(i)
               Call Quit(_RC_CHO_RUN_)
            Else
               Diag(i)=zero
            End If
         End If
      End Do

      Return
      End
