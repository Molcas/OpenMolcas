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
! Copyright (C) 2006, Francesco Aquilante                              *
!               2014, Thomas Bondo Pedersen                            *
!***********************************************************************
!  INV_CHO_FACTOR
!
!> @brief
!>   Evaluation of the inverse Cholesky factor (\f$ Q \f$) of a SPD matrix (\f$ A \f$)
!>   by using a modified Gram--Schmidt orthonormalization of a set
!>   of unit vectors
!> @author F. Aquilante (Nov. 2006)
!> @modified_by T.B. Pedersen (2014) Change criterion for too negative norm
!>
!> @details
!> Evaluation of the inverse Cholesky factor (\f$ Q \f$) of a SPD matrix (\f$ A \f$)
!> by using a modified Gram--Schmidt orthonormalization of a set of unit vectors \f$ V: V(i,k)=\delta_{ik} \f$):
!>
!> \code
!>   For k=1,dim(A)
!>     Qu_k = V_k - sum_j=1^k-1 (Q_j^T * A * V_k) * Q_j
!>          = V_k - sum_j=1^k-1 (Q_j^T * A_k) * Q_j
!>     Q_k = Qu_k / sqrt(Qu_k^T * A_k * Qu_k)
!> \endcode
!>
!> The result is such that the inverse of A is Cholesky decomposed as
!>
!> \f[ A^{-1} = Q Q^\text{T} \f]
!> (\f$ Q \f$ is a full-rank upper triangular matrix)
!>
!> or in general (also for rank deficient \f$ A \f$) such that
!>
!> \f[ Q^\text{T} A Q = I \f]
!>
!> The inverse Cholesky factor is in general *NOT UNIQUE!!*
!> Therefore, and for stability reason, a full pivoting of the
!> initial matrix \f$ A \f$ would be advisable.
!>
!> Worth of mention is the fact that the lower triangular
!> matrix L such that
!>
!> \f[ A = L L^\text{T} \f]
!> (Cholesky decomposition of \f$ A \f$)
!>
!> can be computed as: \f$ L = A Q \f$
!>
!> @side_effects
!> In output \p A_k is returned in a PACKED form (i.e. off-diagonal elements are
!> scaled by two); the latter is the form in which it should be stored as column of \p Am.
!> In case of detected linear dependence, the \p Q_k array is returned as zeros!
!>
!> @note
!> Triangular storage must be used for the \f$ Q \f$-matrix!
!>
!> @param[in,out] A_k    \p kCol -th column of \f$ A \f$ (min. size \p kCol)
!> @param[in]     kCol   index of the column/vector
!> @param[in]     Am     in-core part of the matrix \f$ A \f$ (triangular storage)
!> @param[in]     Qm     in-core matrix whose columns are the orthonormal vectors (triangular storage)
!> @param[in]     nMem   max number of columns of \p Qm (and also of \p Am) kept in core
!> @param[in]     lu_A   file unit where the \f$ A \f$-matrix is stored
!> @param[in]     lu_Q   file unit where the \f$ Q \f$-matrix is stored
!> @param[in]     Scr    scratch space used for reading out-of-core columns of \p Qm and \p Am
!> @param[in]     lScr   size of the scratch space (&ge; \p kCol-1 or ``0`` iff in-core)
!> @param[in]     Z      auxiliary array of min. size \p kCol (always needed)
!> @param[in]     X      auxiliary array of min. size \p kCol-1 (needed only for the out-of-core case)
!> @param[in]     thr    threshold for linear dependence
!> @param[out]    Q_k    the \p kCol -th column of \p Qm (min. size \p kCol)
!> @param[out]    lindep integer indicating detected linear dependence (= ``1`` iff found lin dep, else = ``0``)
!***********************************************************************

subroutine INV_CHO_FACTOR(A_k,kCol,Am,Qm,nMem,lu_A,lu_Q,Scr,lScr,Z,X,thr,Q_k,lindep)

use Index_Functions, only: iTri, nTri_Elem
#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs, Is_Real_Par
#endif
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: A_k(*)
integer(kind=iwp), intent(in) :: kCol, nMem, lu_A, lu_Q, lScr
real(kind=wp), intent(in) :: Am(*), Qm(*), thr
real(kind=wp), intent(_OUT_) :: Scr(*), Z(*), X(*), Q_k(*)
integer(kind=iwp), intent(out) :: lindep
#include "warnings.h"
integer(kind=iwp) :: i, IJ, j, jp, kdone, kread, kstart, lQcol, lQdone, lQdone_, lQread
real(kind=wp) :: sprev, xnorm
real(kind=wp), parameter :: thr_neg = -1.0e-8_wp
real(kind=wp), external :: ddot_

!***********************************************************************
if (thr < zero) then
  call WarningMessage(2,'Error in Inv_Cho_Factor')
  write(u6,*) 'thr must be >= zero'
  call Quit(_RC_CHO_LOG_)
end if

lindep = 0

if (kCol <= nMem) then

  ! Compute scalar product of A_k with previous vectors
  ! ---------------------------------------------------
  jp = 1
  do j=1,kCol-1
    Z(j) = ddot_(j,A_k(1),1,Qm(jp),1)
    jp = jp+j
  end do
  !call RecPrt('A_k*Qm',' ',Z,1,kCol)

  ! Compute unnormalized k-th vector
  ! --------------------------------
  ! SVC: this piece of code was computing Q_k = - Qm * Z, where Q_k and Z
  ! are vectors of length kCol-1, and Qm is a matrix in triangular storage
  ! with column-wise layout:
  !   |Q_k(1)     |      |Qm(1,1) Qm(1,2) ... Qm(1,kCol-1)     |   |Z(1)     |
  !   |Q_k(2)     |      |        Qm(2,2) ... Qm(2,kCol-1)     |   |Z(2)     |
  !   | ...       |  = - |                ...     ...          | * |...      |
  !   | ...       |      |                                     |   |...      |
  !   |Q_k(kCol-1)|      |                    Qm(kCol-1,kCol-1)|   |Z(kCol-1)|
  ! In order to improve performance, I've used the DTPMV routine from
  ! BLAS. For parallel processes, we will block up the triangular matrix
  ! and divide the blocks over the processes, using either DTPMV or DGEMV
  ! on the blocks (depending if it is a diagonal or off-diagonal block).

  Q_k(1:kCol-1) = Zero
# ifdef _MOLCAS_MPP_
  if (is_real_par() .and. (kCol >= 500)) then
    ! SVC: the best way would probably be to chop up the triangular matrix
    ! into blocks, and then call DTPMV/DGEMV on those blocks. To keep things
    ! simple, I've just used a series of DAXPY's on each column of the
    ! triangular matrix, this should be sufficient (for now).
    do J=1+MYRANK,KCOL-1,NPROCS
      IJ = iTri(J,1)
      Q_k(1:J) = Q_k(1:J)-Z(J)*Qm(IJ:IJ+J-1)
    end do
    call GAdGOp(Q_k,kCol-1,'+')
  else
# endif
    Q_k(1:kCol-1) = Q_k(1:kCol-1)-Z(1:kCol-1)
    call DTPMV('U','N','N',kCol-1,Qm,Q_k,1)
# ifdef _MOLCAS_MPP_
  end if
# endif
  Q_k(kCol) = one

  ! Normalize k-th vector :   ||Q_k|| = Q_k^T * A * Q_k
  ! ---------------------------------------------------

  A_k(1:kCol-1) = Two*A_k(1:kCol-1) ! packing of A_k

  Z(kCol) = ddot_(kCol,A_k(1),1,Q_k(1),1) !contrib fr k-th col of A

  jp = 1
  do j=1,kCol-1 ! contrib. from previous columns of A
    Z(j) = ddot_(j,Q_k(1),1,Am(jp),1)
    jp = jp+j
  end do

  xnorm = ddot_(kCol,Z(1),1,Q_k(1),1)

  if (xnorm >= thr) then

    xnorm = one/sqrt(xnorm)
    Q_k(1:kCol) = xnorm*Q_k(1:kCol)

    !-tbp: use fixed criterion for too negative diagonal
    !else if ((xnorm > zero) .or. (-xnorm <= Ten*thr)) then
  else if (xnorm > thr_neg) then

    lindep = 1
    Q_k(1:kCol) = Zero

  else

    call WarningMessage(2,'Error in Inv_Cho_Factor')
    write(u6,*) 'INV_CHO_FACTOR: too-negative value for norm(Q_k).'
    write(u6,*) 'INV_CHO_FACTOR: xnorm = ',xnorm
    call Quit(_RC_CHO_RUN_)

  end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
else   ! the first nMem columns of Q are in memory
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  if (lScr < kCol-1) then
    call WarningMessage(2,'Error in Inv_Cho_Factor')
    write(u6,*) 'lScr must be >= kCol-1'
    call Quit(_RC_CHO_LOG_)
  end if

  X(1:kCol-1) = Zero

  ! Compute scalar product of A_k with in-core previous vectors
  ! -----------------------------------------------------------
  jp = 1
  do j=1,nMem
    Z(j) = ddot_(j,A_k(1),1,Qm(jp),1)
    jp = jp+j
  end do

  ! Batch for the out-of-core previous vectors
  ! ------------------------------------------
  kdone = nMem
  lQcol = nTri_Elem(kCol-1) ! length up to kCol-1
  do while (kdone < kCol-1)

    lQdone = nTri_Elem(kdone)
    lQdone_ = lQdone
    lQread = lQcol-lQdone

    kread = kCol-1
    do while (lQread > lScr)
      lQread = lQread-kread
      kread = kread-1
    end do

    call ddafile(lu_Q,2,Scr,lQread,lQdone_) ! read

    jp = 1
    do j=kdone+1,kread
      Z(j) = ddot_(j,A_k,1,Scr(jp),1)
      jp = jp+j
    end do

    ! Store an out-of-core intermediate for the Q-vectors
    ! ---------------------------------------------------
    do i=1,kread
      sprev = zero
      kstart = max(i,kdone+1) ! ((j >= i) .and. j_out_of_core)
      do j=kstart,kread
        ij = iTri(j,i)-lQdone
        sprev = sprev+Z(j)*Scr(ij)
      end do
      X(i) = X(i)+sprev
    end do

    kdone = kread

  end do
  !call RecPrt('A_k*Qm',' ',Z,1,kCol)

  ! Compute unnormalized k-th vector
  ! --------------------------------
  do i=1,kCol-1
    sprev = X(i) ! out-of-core contrib.
    do j=i,nMem
      ij = iTri(j,i)
      sprev = sprev+Z(j)*Qm(ij)
    end do
    Q_k(i) = -sprev
  end do
  Q_k(kCol) = one

  ! Normalize k-th vector :   ||Q_k|| = Q_k^T * A * Q_k
  ! ---------------------------------------------------

  A_k(1:kCol-1) = Two*A_k(1:kCol-1) ! packing of A_k

  Z(kCol) = ddot_(kCol,A_k(1),1,Q_k(1),1) !contrib fr k-th col of A

  ! Batch for the out-of-core previous vectors
  ! ------------------------------------------
  kdone = nMem
  lQcol = nTri_Elem(kCol-1) ! length up to kCol-1
  do while (kdone < kCol-1)

    lQdone = nTri_Elem(kdone)
    lQread = lQcol-lQdone

    kread = kCol-1
    do while (lQread > lScr)
      lQread = lQread-kread
      kread = kread-1
    end do

    ! Out-of-core intermediate to be used for the normalization factor
    ! ----------------------------------------------------------------
    call ddafile(lu_A,2,Scr,lQread,lQdone) ! read

    jp = 1
    do j=kdone+1,kread
      Z(j) = ddot_(j,Q_k,1,Scr(jp),1)
      jp = jp+j
    end do

    kdone = kread

  end do

  jp = 1
  do j=1,nMem ! contrib. from in-core previous columns of A
    Z(j) = ddot_(j,Q_k(1),1,Am(jp),1)
    jp = jp+j
  end do

  xnorm = ddot_(kCol,Z(1),1,Q_k(1),1)

  if (xnorm >= thr) then

    xnorm = one/sqrt(xnorm)
    Q_k(1:kCol) = xnorm*Q_k(1:kCol)

    !-tbp: use fixed criterion for too negative diagonal
    !else if ((xnorm > zero) .or. (-xnorm <= Ten*thr)) then
  else if (xnorm > thr_neg) then

    lindep = 1
    Q_k(1:kCol) = Zero

  else

    call WarningMessage(2,'Error in Inv_Cho_Factor')
    write(u6,*) 'INV_CHO_FACTOR: too-negative value for norm(Q_k).'
    write(u6,*) 'INV_CHO_FACTOR: xnorm = ',xnorm
    call Quit(_RC_CHO_RUN_)

  end if

end if

return

end subroutine INV_CHO_FACTOR
