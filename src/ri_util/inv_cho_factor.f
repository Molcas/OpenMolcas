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
* Copyright (C) 2006, Francesco Aquilante                              *
*               2014, Thomas Bondo Pedersen                            *
************************************************************************
*  INV_CHO_FACTOR
*
*> @brief
*>   Evaluation of the inverse Cholesky factor (\f$ Q \f$) of a SPD matrix (\f$ A \f$)
*>   by using a modified Gram--Schmidt orthonormalization of a set
*>   of unit vectors
*> @author F. Aquilante (Nov. 2006)
*> @modified_by T.B. Pedersen (2014) Change criterion for too negative norm
*>
*> @details
*> Evaluation of the inverse Cholesky factor (\f$ Q \f$) of a SPD matrix (\f$ A \f$)
*> by using a modified Gram--Schmidt orthonormalization of a set of unit vectors \f$ V: V(i,k)=\delta_{ik} \f$):
*>
*> \code
*>   For k=1,dim(A)
*>     Qu_k = V_k - sum_j=1^k-1 (Q_j^T * A * V_k) * Q_j
*>          = V_k - sum_j=1^k-1 (Q_j^T * A_k) * Q_j
*>     Q_k = Qu_k / sqrt(Qu_k^T * A_k * Qu_k)
*> \endcode
*>
*> The result is such that the inverse of A is Cholesky decomposed as
*>
*> \f[ A^{-1} = Q Q^\text{T} \f]
*> (\f$ Q \f$ is a full-rank upper triangular matrix)
*>
*> or in general (also for rank deficient \f$ A \f$) such that
*>
*> \f[ Q^\text{T} A Q = I \f]
*>
*> The inverse Cholesky factor is in general *NOT UNIQUE!!*
*> Therefore, and for stability reason, a full pivoting of the
*> initial matrix \f$ A \f$ would be advisable.
*>
*> Worth of mention is the fact that the lower triangular
*> matrix L such that
*>
*> \f[ A = L L^\text{T} \f]
*> (Cholesky decomposition of \f$ A \f$)
*>
*> can be computed as: \f$ L = A Q \f$
*>
*> @side_effects
*> In output \p A_k is returned in a PACKED form (i.e. off-diagonal elements are
*> scaled by two); the latter is the form in which it should be stored as column of \p Am.
*> In case of detected linear dependence, the \p Q_k array is returned as zeros!
*>
*> @note
*> Triangular storage must be used for the \f$ Q \f$-matrix!
*>
*> @param[in,out] A_k    \p kCol -th column of \f$ A \f$ (min. size \p kCol)
*> @param[in]     kCol   index of the column/vector
*> @param[in]     Am     in-core part of the matrix \f$ A \f$ (triangular storage)
*> @param[in]     Qm     in-core matrix whose columns are the orthonormal vectors (triangular storage)
*> @param[in]     nMem   max number of columns of \p Qm (and also of \p Am) kept in core
*> @param[in]     lu_A   file unit where the \f$ A \f$-matrix is stored
*> @param[in]     lu_Q   file unit where the \f$ Q \f$-matrix is stored
*> @param[in]     Scr    scratch space used for reading out-of-core columns of \p Qm and \p Am
*> @param[in]     lScr   size of the scratch space (&ge; \p kCol-1 or ``0`` iff in-core)
*> @param[in]     Z      auxiliary array of min. size \p kCol (always needed)
*> @param[in]     X      auxiliary array of min. size \p kCol-1 (needed only for the out-of-core case)
*> @param[in]     thr    threshold for linear dependence
*> @param[out]    Q_k    the \p kCol -th column of \p Qm (min. size \p kCol)
*> @param[out]    lindep integer indicating detected linear dependence (= ``1`` iff found lin dep, else = ``0``)
************************************************************************
      SUBROUTINE INV_CHO_FACTOR(A_k,kCol,Am,Qm,nMem,lu_A,lu_Q,Scr,lScr,
     &                          Z,X,thr,Q_k,lindep)

      Implicit Real*8 (a-h,o-z)
      Integer kCol, nMem, lu_A, lu_Q, lScr, lindep
      Integer  incr,shft
      Real*8  A_k(*), Am(*), Qm(*), Scr(*), Z(*), X(*), Q_k(*)
      Real*8  thr
#include "para_info.fh"
#include "warnings.fh"
      Parameter ( two = 2.0d0, one = 1.0d0, zero = 0.0d0 )
      Parameter ( thr_neg = -1.0d-8 )
**********************************************************************
      incr=1
      shft=0

      If (thr .lt. zero) Then
         Call WarningMessage(2,'Error in Inv_Cho_Factor')
         write(6,*)'thr must be .ge. zero'
         Call Quit(_RC_CHO_LOG_)
      EndIf

      lindep = 0

      If (kCol .le. nMem) Then
*
*  Compute scalar product of A_k with previous vectors
*  ---------------------------------------------------
         jp=1
         Do j=1,kCol-1
            Z(j)=ddot_(j,A_k(1),1,Qm(jp),1)
            jp = jp + j
         End Do
C        Call RecPrt('A_k*Qm',' ',Z,1,kCol)
*
*  Compute unnormalized k-th vector
* ---------------------------------
C SVC: this piece of code was computing Q_k = - Qm * Z, where Q_k and Z
C are vectors of length kCol-1, and Qm is a matrix in triangular storage
C with column-wise layout:
C   |Q_k(1)     |      |Qm(1,1) Qm(1,2) ... Qm(1,kCol-1)     |   |Z(1)     |
C   |Q_k(2)     |      |        Qm(2,2) ... Qm(2,kCol-1)     |   |Z(2)     |
C   | ...       |  = - |                ...     ...          | * |...      |
C   | ...       |      |                                     |   |...      |
C   |Q_k(kCol-1)|      |                    Qm(kCol-1,kCol-1)|   |Z(kCol-1)|
C In order to improve performance, I've used the DTPMV routine from
C BLAS. For parallel processes, we will block up the triangular matrix
C and divide the blocks over the processes, using either DTPMV or DGEMV
C on the blocks (depending if it is a diagonal or off-diagonal block).

         Call FZero(Q_k,kCol-1)
#ifdef _MOLCAS_MPP_
         if (is_real_par().and.kCol.ge.500) then
C SVC: the best way would probably be to chop up the triangular matrix
C into blocks, and then call DTPMV/DGEMV on those blocks. To keep things
C simple, I've just used a series of DAXPY's on each column of the
C triangular matrix, this should be sufficient (for now).
           DO J=1+MYRANK,KCOL-1,NPROCS
             IJ=(J*(J-1))/2+1
             CALL DAXPY_(J,-Z(J),Qm(IJ),1,Q_k,1)
           END DO
           Call GAdGOp(Q_k,kCol-1,'+')
         else
           CALL DAXPY_(kCol-1,-1.0D0,Z,1,Q_k,1)
           CALL DTPMV('U','N','N',kCol-1,Qm,Q_k,1)
         end if
#else
         CALL DAXPY_(kCol-1,-1.0D0,Z,1,Q_k,1)
         CALL DTPMV('U','N','N',kCol-1,Qm,Q_k,1)
#endif
         Q_k(kCol) = one
*
*
*  Normalize k-th vector :   ||Q_k|| = Q_k^T * A * Q_k
* ----------------------------------------------------

         call dscal_(kCol-1,two,A_k(1),1) ! packing of A_k

         Z(kCol)=ddot_(kCol,A_k(1),1,Q_k(1),1) !contrib fr k-th col of A

         jp=1
         Do j=1,kCol-1 ! contrib. from previous columns of A
            Z(j)=ddot_(j,Q_k(1),1,Am(jp),1)
            jp = jp + j
         End Do

         xnorm=ddot_(kCol,Z(1),1,Q_k(1),1)

         If (xnorm.ge.thr) Then

           xnorm=one/sqrt(xnorm)
           call dscal_(kCol,xnorm,Q_k(1),1)

C-tbp: use fixed criterion for too negative diagonal
C-tbp    ElseIf (xnorm.gt.zero .or. -xnorm.le.1.0d1*thr) Then
         ElseIf (xnorm.gt.thr_neg) Then

           lindep = 1
           Call Fzero(Q_k(1),kCol)

         Else

           Call WarningMessage(2,'Error in Inv_Cho_Factor')
           write(6,*)'INV_CHO_FACTOR: too-negative value for norm(Q_k).'
           write(6,*)'INV_CHO_FACTOR: xnorm = ',xnorm
           Call Quit(_RC_CHO_RUN_)

         EndIf


*                                                                      *
************************************************************************
*                                                                      *
      Else   ! the first nMem columns of Q are in memory
*                                                                      *
************************************************************************
*                                                                      *

         If (lScr .lt. kCol-1) Then
            Call WarningMessage(2,'Error in Inv_Cho_Factor')
            write(6,*)'lScr must be .ge. kCol-1'
            Call Quit(_RC_CHO_LOG_)
         EndIf

         Call FZero(X(1),kCol-1)
*
*  Compute scalar product of A_k with in-core previous vectors
*  -----------------------------------------------------------
         jp=1
         Do j=1,nMem
            Z(j)=ddot_(j,A_k(1),1,Qm(jp),1)
            jp = jp + j
         End Do
*
*  Batch for the out-of-core previous vectors
*--------------------------------------------
         kdone = nMem
         lQcol = (kCol-1)*kCol/2 ! length up to kCol-1
         Do while ( kdone .lt. kCol-1 )

            lQdone = kdone*(kdone+1)/2
            lQdone_=lQdone
            lQread = lQcol - lQdone

            kread = kCol-1
            Do while ( lQread .gt. lScr )
              lQread = lQread - kread
              kread = kread - 1
            End Do

            Call ddafile(lu_Q,2,Scr(1),lQread,lQdone_) ! read

            jp=1
            Do j=kdone+1,kread
               Z(j)=ddot_(j,A_k(1),1,Scr(jp),1)
               jp = jp + j
            End Do
*
*  Store an out-of-core intermediate for the Q-vectors
*-----------------------------------------------------
            Do i=1,kread
              sprev=zero
              kstart = Max(i,kdone+1) ! (j.ge.i .and. j_out_of_core)
              Do j=kstart,kread
                 ij = j*(j-1)/2 + i - lQdone
                 sprev = sprev + Z(j)*Scr(ij)
              End Do
              X(i) = X(i) + sprev
            End Do

            kdone = kread

         End Do
C        Call RecPrt('A_k*Qm',' ',Z,1,kCol)
*
*  Compute unnormalized k-th vector
* ---------------------------------
         Do i=1,kCol-1
            sprev = X(i) ! out-of-core contrib.
            Do j=i,nMem
               ij = j*(j-1)/2 + i
               sprev = sprev + Z(j)*Qm(ij)
            End Do
            Q_k(i) = - sprev
         End Do
         Q_k(kCol) = one
*
*
*  Normalize k-th vector :   ||Q_k|| = Q_k^T * A * Q_k
* ----------------------------------------------------

         call dscal_(kCol-1,two,A_k(1),1) ! packing of A_k

         Z(kCol)=ddot_(kCol,A_k(1),1,Q_k(1),1) !contrib fr k-th col of A
*
*  Batch for the out-of-core previous vectors
*--------------------------------------------
         kdone = nMem
         lQcol = (kCol-1)*kCol/2 ! length up to kCol-1
         Do while ( kdone .lt. kCol-1 )

            lQdone = kdone*(kdone+1)/2
            lQread = lQcol - lQdone

            kread = kCol-1
            Do while ( lQread .gt. lScr )
              lQread = lQread - kread
              kread = kread - 1
            End Do
*
*  Out-of-core intermediate to be used for the normalization factor
*------------------------------------------------------------------
            Call ddafile(lu_A,2,Scr(1),lQread,lQdone) ! read

            jp=1
            Do j=kdone+1,kread
               Z(j)=ddot_(j,Q_k(1),1,Scr(jp),1)
               jp = jp + j
            End Do

            kdone = kread

         End Do

         jp=1
         Do j=1,nMem ! contrib. from in-core previous columns of A
            Z(j)=ddot_(j,Q_k(1),1,Am(jp),1)
            jp = jp + j
         End Do

         xnorm=ddot_(kCol,Z(1),1,Q_k(1),1)

         If (xnorm.ge.thr) Then

           xnorm=one/sqrt(xnorm)
           call dscal_(kCol,xnorm,Q_k(1),1)

C-tbp: use fixed criterion for too negative diagonal
C-tbp    ElseIf (xnorm.gt.zero .or. -xnorm.le.1.0d1*thr) Then
         ElseIf (xnorm.gt.thr_neg) Then

           lindep = 1
           Call Fzero(Q_k(1),kCol)

         Else

           Call WarningMessage(2,'Error in Inv_Cho_Factor')
           write(6,*)'INV_CHO_FACTOR: too-negative value for norm(Q_k).'
           write(6,*)'INV_CHO_FACTOR: xnorm = ',xnorm
           Call Quit(_RC_CHO_RUN_)

         EndIf


      EndIf

      Return
      End
