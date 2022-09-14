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
! Copyright (C) 2007, Francesco Aquilante                              *
!               2014, Thomas Bondo Pedersen                            *
!***********************************************************************
!  CHO_FACTOR
!
!> @brief
!>   Evaluation of the Cholesky factor (\f$ Z \f$) of a SPD matrix (\f$ A \f$)
!> @author F. Aquilante (Jan. 2007)
!> @modified_by T.B. Pedersen (2014) Changed criterion for too negative diagonal
!>
!> @details
!> Evaluation of the Cholesky factor (\f$ Z \f$) of a SPD matrix (\f$ A \f$)
!>
!> \code
!>   For k=1,dim(A)
!>     Z(k,k) = sqrt( A(k,k) - sum_j  Z(k,j)^2 )
!>     Z(i,k) = ( A(i,k) - sum_j  Z(i,j)*Z(k,j) ) / Z(k,k)
!> \endcode
!>
!> The result is such that \f$ A \f$ is Cholesky decomposed as
!>
!> \f[  A = Z Z^\text{T} \f]
!>
!> The Cholesky factor is in general *NOT UNIQUE!!*
!> Therefore, and for stability reason, pivoting of the
!> initial matrix \f$ A \f$ would be advisable.
!>
!> @side_effects
!> \p A_k in output contains the \p kCol -th Cholesky vector.
!> In case of detected linear dependence, the \p A_k array
!> is returned as zeros!
!>
!> @note
!> Rectangular storage must be used for the \f$ Z \f$-matrix!
!>
!> @param[in,out] Diag   Updated diagonal elements of \f$ A \f$ (subtraction done by this routine)
!> @param[in,out] A_k    currently treated column of \f$ A \f$. In output contains the \p kCol -th Cholesky vector
!> @param[in]     iD_A   indices of the columns of \f$ A \f$
!> @param[in]     kCol   index of the Cholesky vector
!> @param[in]     nRow   number of rows of \f$ A \f$
!> @param[in]     Zm     in-core matrix whose columns are the Cholesky vectors
!> @param[in]     nMem   max number of columns of \p Zm kept in core
!> @param[in]     lu_Z   file unit where the \f$ Z \f$-matrix is stored
!> @param[out]    Scr    scratch space used for reading out-of-core columns of \p Zm
!> @param[in]     lScr   size of the scratch space (&ge; \p nRow or ``0`` iff in-core)
!> @param[in]     thr    threshold for linear dependence
!> @param[out]    lindep integer indicating detected linear dependence (= ``1`` iff found lin dep, else = ``0``)
!***********************************************************************

subroutine CHO_FACTOR(Diag,A_k,iD_A,kCol,nRow,Zm,nMem,lu_Z,Scr,lScr,thr,lindep)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: Diag(*), A_k(*)
integer(kind=iwp), intent(in) :: iD_A(*), kCol, nRow, nMem, lu_Z, lScr
real(kind=wp), intent(in) :: Zm(nRow,*), thr
real(kind=wp), intent(_OUT_) :: Scr(*)
integer(kind=iwp), intent(out) :: lindep
#include "warnings.h"
integer(kind=iwp) :: i, j, kdone, kj, kstep, lZdone, lZread, lZrem
real(kind=wp) :: Dmax, fac, xfac
real(kind=wp), parameter :: thr_neg = -1.0e-8_wp

!***********************************************************************

if (thr < zero) then
  call WarningMessage(2,'Error in Cho_Factor')
  write(u6,*) 'thr must be >= zero'
  call Quit(_RC_CHO_LOG_)
end if

lindep = 0
Dmax = Diag(iD_A(kCol)) ! pivoting done by the calling routine
xfac = one/sqrt(abs(Dmax))

if (kCol <= nMem) then

  if (Dmax >= thr) then

    ! Compute elements of the k-th Cholesky vector
    !---------------------------------------------
    !    Z(i,k) = A(i,k) - sum_j  Z(k,j)*Z(i,j)
    !---------------------------------------------
    do j=1,kCol-1

      fac = -Zm(iD_A(kCol),j)
      A_k(1:nRow) = A_k(1:nRow)+fac*Zm(:,j)

    end do

    !-tbp: use thr_neg as threshold for too negative diagonal
    !      It should not depend on the decomposition threshold!
  !else if ((Dmax > Zero) .or. (-Dmax <= Ten*thr)) then
  else if (Dmax > thr_neg) then

    lindep = 1
    A_k(1:nRow) = Zero
    return

  else

    call WarningMessage(2,'Error in Cho_Factor')
    write(u6,*) 'CHO_FACTOR: too-negative diagonal.'
    write(u6,*) 'CHO_FACTOR: current largest Diag = ',Dmax
    call Quit(_RC_CHO_RUN_)

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else  ! the first nMem columns of Z are in memory
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  if (lScr < nRow) then
    call WarningMessage(2,'Error in Cho_Factor')
    write(u6,*) 'lScr must be >= nRow'
    call Quit(_RC_CHO_LOG_)
  end if

  if (Dmax >= thr) then

    ! Compute elements of the k-th Cholesky vector (in-core contrib.)
    !----------------------------------------------------------------
    !    Z(i,k) = A(i,k) - sum_j  Z(k,j)*Z(i,j)
    !----------------------------------------------------------------
    do j=1,nMem

      fac = -Zm(iD_A(kCol),j)
      A_k(1:nRow) = A_k(1:nRow)+fac*Zm(:,j)

    end do

    ! Batch for the out-of-core previous vectors
    !-------------------------------------------
    kstep = lScr/nRow

    do kdone=nMem+1,kCol-1,kStep

      lZdone = nRow*(kdone-1)
      lZrem = nRow*(kCol-kdone)
      lZread = min(LZrem,nRow*kStep)

      call ddafile(lu_Z,2,Scr,lZread,lZdone) ! read

      ! Compute elements of the k-th Cholesky vector (out-of-core contrib.)
      !--------------------------------------------------------------------
      !    Z(i,k) = A(i,k) - sum_j  Z(k,j)*Z(i,j)
      !--------------------------------------------------------------------
      do j=1,lZread/nRow
        kj = nRow*(j-1)
        fac = -Scr(kj+iD_A(kCol))
        A_k(1:nRow) = A_k(1:nRow)+fac*Scr(kj+1:kj+nRow)
      end do

    end do

    !-tbp: use thr_neg as threshold for too negative diagonal
    !      It should not depend on the decomposition threshold!
  !else if ((Dmax > Zero) .or. (-Dmax <= Ten*thr)) then
  else if (Dmax > thr_neg) then

    lindep = 1
    A_k(1:nRow) = Zero
    return

  else

    call WarningMessage(2,'Error in Cho_Factor')
    write(u6,*) 'CHO_FACTOR: too-negative diagonal.'
    write(u6,*) 'CHO_FACTOR: current largest Diag = ',Dmax
    call Quit(_RC_CHO_RUN_)

  end if

end if

A_k(iD_A(kCol)) = Dmax

! Scaling of the vector elements :  Z(i,k) = Z(i,k)/Z(k,k)
! --------------------------------------------------------
A_k(1:nRow) = xfac*A_k(1:nRow)

!  Explicit zeroing of the previously treated elements
! ----------------------------------------------------
do i=1,kCol-1
  A_k(iD_A(i)) = zero
end do

! Update diagonal elements of the A matrix
!------------------------------------------
!   A(i,i) = A(i,i) - Z(i,k)^2    ( i > k )
!------------------------------------------
Diag(1:nRow) = Diag(1:nRow)-A_k(1:nRow)**2
Diag(iD_A(kCol)) = zero ! explicit zeroing of the treated diagonal

!-tbp: zero negative diagonal elements
!      Stop if too negative!
do i=1,nRow
  if (Diag(i) < zero) then
    if (Diag(i) <= thr_neg) then
      call WarningMessage(2,'Error in Cho_Factor')
      write(u6,*) 'CHO_FACTOR: too negative diagonal.'
      write(u6,*) 'CHO_FACTOR: i,Diag(i)= ',i,Diag(i)
      call Quit(_RC_CHO_RUN_)
    else
      Diag(i) = zero
    end if
  end if
end do

return

end subroutine CHO_FACTOR
