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
! Copyright (C) Thomas Bondo Pedersen                                  *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  LinEqSolv
!
!> @brief
!>   Solve linear equations \f$ Ax=B \f$ or \f$ A^\text{T}x=B \f$ where
!>   \f$ A \f$ is a general nonsingular matrix
!> @author Thomas Bondo Pedersen
!> @modified_by Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> For \p TransA(1:1) = ``'N'`` or ``'n'``,
!> this routine solves the equations \f$ Ax = B \f$ for any number of
!> righthand sides stored as columns of the matrix \f$ B \f$. For
!> \p TransA(1:1) = ``'T'`` or ``'t'``, the equations \f$ A^\text{T}x = B \f$ are
!> solved.
!>
!> The solution vectors are stored as columns of \p B on exit. Note
!> that array \p A will be destroyed during the solution (the matrix
!> is replaced by its factors obtained by Gaussian elimination).
!>
!> Return codes are:
!>
!> - \p irc = ``-1``: input error(s) detected and nothing has been done.
!> - \p irc =  ``0``: successful completion.
!> - \p irc =  ``1``: \p A is estimated to be singular and no solution vectors have been computed.
!>
!> @note
!> 4 * \p nDim real*8 and 2 * \p nDim integer words of memory must be available
!> on entry.
!>
!> @param[out]    irc    Return code
!> @param[in]     TransA Transposition of \p A
!> @param[in,out] A      Array containing the nonsingular matrix \p A on entry and the factorized matrix (Gaussian elimination)
!>                       on exit
!> @param[in]     ldA    Leading dimension of \p A
!> @param[in,out] B      Array containing righthand sides of equations on entry and solution vectors on exit
!> @param[in]     ldB    Leading dimension of \p B
!> @param[in]     nDim   Column dimension of \p A
!> @param[in]     nEq    Number of equations, i.e. column dimension of \p B
!***********************************************************************

subroutine LinEqSolv(irc,TransA,A,ldA,B,ldB,nDim,nEq)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
character(len=*), intent(in) :: TransA
integer(kind=iwp), intent(in) :: ldA, ldB, nDim, nEq
real(kind=wp), intent(inout) :: A(ldA,nDim), B(ldB,nEq)
integer(kind=iwp) :: iErr, lTransA
real(kind=wp) :: AN, RC
character :: myTransA
integer(kind=iwp), allocatable :: iPivot(:), iScr(:)
real(kind=wp), allocatable :: Scr(:)
real(kind=wp), external :: DLANGE_

! Test input.
! -----------

irc = 0
if ((nEq < 1) .or. (nDim < 1)) return
if ((ldA < nDim) .or. (ldB < nDim)) then
  irc = -1 ! dimension error
  return
end if

! Translate TransA to integer Job to be used by DGESL.
! ----------------------------------------------------

lTransA = len(TransA)
if (lTransA > 0) then
  myTransA = TransA(1:1)
else
  irc = -1 ! TransA error
  return
end if

if ((myTransA /= 'N') .and. (myTransA /= 'n') .and. (myTransA /= 'T') .and. (myTransA /= 't')) then
  irc = -1 ! TransA error
  return
end if

! Allocate pivot array and scratch space
! --------------------------------------

call mma_allocate(iPivot,nDim,label='LES_Pivot')
call mma_allocate(Scr,4*nDim,label='LES_Scr')
call mma_allocate(iScr,nDim,label='LES_iScr')

! Factor A by Gaussian elimination and estimate reciprocal condition
! number (RC). Check for singularity.
! ------------------------------------------------------------------

RC = Zero
AN = DLANGE_('1',nDim,nDim,A,ldA,Scr)
call DGETRF_(nDim,nDim,A,ldA,iPivot,iErr)
call DGECON_('1',nDim,A,ldA,AN,RC,Scr,iScr,iErr)
if ((One+RC == One) .or. (iErr > 0)) then
  irc = 1 ! error: A is (probably) singular
else

  ! Solve equations.
  ! ----------------

  call DGETRS_(myTransA,nDim,nEq,A,ldA,iPivot,B,ldB,iErr)
  if (iErr > 0) irc = 1 ! error: A is (probably) singular

end if

! Deallocations.
! --------------

call mma_deallocate(iPivot)
call mma_deallocate(Scr)
call mma_deallocate(iScr)

end subroutine LinEqSolv
