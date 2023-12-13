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
!***********************************************************************
!  ChoDec
!
!> @brief
!>   Decompose a symmetric positive (semi-)definite matrix
!> @author Thomas Bondo Pedersen
!>
!> @details
!> ::ChoDec will Cholesky decompose a positive (semi-)definite
!> matrix without knowledge of all its elements (i.e., using a
!> matrix-direct approach).
!>
!> Given the diagonal, the idea is that ::ChoDec asks the caller
!> for specified columns of the matrix when needed during
!> decomposition. This is done by invoking the external routine
!> passed as \p CD_Col. Similarly, the caller must provide an
!> external routine, passed as \p CD_Vec, that handles the Cholesky
!> vector buffer---that is, empties or fills the buffer (\p CD_Vec
!> could e.g. do Cholesky vector I/O).
!>
!> The external subroutine \p CD_Col must have the following call
!> pattern:
!>
!> \code
!>   Call CD_Col(Col,nDim,iCol,nCol,Buf,lBuf)
!> \endcode
!>
!> where \p Col(nDim,nCol) contains the columns, \p iCol(nCol) specifies
!> the columns requested [i.e., column \p Col(*,i) should correspond
!> to column \p iCol(i) in the full matrix], and \p Buf(lBuf) can be
!> used as (perhaps additional) work space.
!>
!> The external subroutine \p CD_Vec must have the following call
!> pattern:
!>
!> \code
!>   Call CD_Vec(iVec1,nVec,Buf,lBuf,nDim,iOpt)
!> \endcode
!>
!> where \p iVec1 is the index of the first vector, \p nVec the number
!> of vectors, and \p nDim the vector dimension. The buffer \p Buf has
!> total dimension \p lBuf. For \p iOpt = ``1`` ("empty buffer option"),
!> the first \p nDim*\p nVec elements of \p Buf contain the vectors on
!> input---i.e. \p CD_Vec is expected to move the vectors to another
!> location so that all of \p Buf is again available for use.
!> For \p iOpt = ``2`` ("fill buffer option"), the process is reversed and
!> \p ChoDec expects \p CD_Vec to return vectors
!> \p iVec1, \p iVec1+1, ..., \p iVec1+\p nVec-1, in the first \p nDim*\p nVec elements
!> of \p Buf.
!>
!> The error code \p irc is ``0`` if success, while a non-zero number
!> is returned if some failure (or if the matrix is not positive
!> (semi-)definite) occurs. In the latter cases, the
!> decomposition is ill-defined and the data returned is junk!
!>
!> @note
!> \p Thr and \p Span may be reset in this routine (if they
!> are not properly defined on input).
!>
!> @param[in]     CD_Col  External subroutine for retrieving matrix columns
!> @param[in]     CD_Vec  External subroutine for filling and emptying Cholesky vector buffer
!> @param[in]     Restart Flag for restarting decomposition
!> @param[in,out] Thr     Decomposition threshold (precision)
!> @param[in,out] Span    Span factor
!> @param[in]     MxQual  Max. number of qualified columns per decomposition pass
!> @param[in]     Diag    Array containing the diagonal of the matrix to be decomposed
!> @param[out]    Qual    Storage array for matrix columns
!> @param[in,out] Buf     Buffer for storing Cholesky vectors
!> @param[out]    iPivot  Pivoting array
!> @param[out]    iQual   Array for storing qualified indices
!> @param[in]     nDim    Linear dimension of matrix
!> @param[in]     lBuf    Buffer size
!> @param[out]    ErrStat Min, max, and RMS errors, respectively, for decomposition as judged
!>                        from the residual diagonal elements
!> @param[in,out] NumCho  Number of Cholesky vectors
!> @param[out]    irc     Return code
!***********************************************************************

subroutine ChoDec(CD_Col,CD_Vec,Restart,Thr,Span,MxQual,Diag,Qual,Buf,iPivot,iQual,nDim,lBuf,ErrStat,NumCho,irc)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
external :: CD_Col, CD_Vec
logical(kind=iwp), intent(in) :: Restart
integer(kind=iwp), intent(in) :: MxQual, nDim, lBuf
real(kind=wp), intent(inout) :: Thr, Span, Buf(lBuf)
real(kind=wp), intent(in) :: Diag(nDim)
real(kind=wp), intent(out) :: Qual(nDim,0:MxQual), ErrStat(3)
integer(kind=iwp), intent(out) :: iPivot(nDim), iQual(MxQual), irc
integer(kind=iwp), intent(inout) :: NumCho
integer(kind=iwp) :: MinBuf, mQual, MxNumCho
logical(kind=iwp) :: Converged
real(kind=wp), parameter :: DefThr = 1.0e-6_wp, DefSpan = 1.0e-2_wp, dum = 9.876543210e15_wp, ThrFail = -1.0e-8_wp, &
                            ThrNeg = -1.0e-13_wp

! Initialize variables.
! ---------------------

irc = 0
ErrStat(1) = dum
ErrStat(2) = -dum
ErrStat(3) = -dum
if (.not. Restart) NumCho = 0
Converged = .false.

! Check dimensions.
! -----------------

if (nDim < 1) then
  irc = 0
  return  ! exit (nothing to do)
end if
if (MxQual < 1) then
  irc = -1
  return  ! exit (qualification not possible)
end if
mQual = min(MxQual,nDim)
MinBuf = nDim+mQual
if (lBuf < MinBuf) then
  irc = -2
  return ! exit (buffer too small)
end if

! Check configuration.
! --------------------

if (Thr < Zero) Thr = DefThr  ! reset threshold
if ((Span < Zero) .or. (Span > One)) Span = DefSpan ! reset span factor

! Set up a (possibly updated) copy of the diagonal.
! Test diagonal for negative elements.
! -------------------------------------------------

call CD_Diag(CD_Vec,Restart,Converged,Thr,ThrNeg,ThrFail,Diag,Qual(:,0),Buf,nDim,lBuf,ErrStat,NumCho,irc)
if (irc /= 0) return ! exit (initial diagonal error)

! Only do decomposition if not converged!
! ---------------------------------------

if (.not. Converged) then

  ! Cholesky decompose.
  ! -------------------

  MxNumCho = nDim ! generate as many vectors as needed
  call CD_Decomposer(CD_Col,CD_Vec,MxNumCho,Thr,Span,mQual,ThrNeg,ThrFail,Qual(:,0),Qual(:,1),Buf,iPivot,iQual,nDim,lBuf,NumCho,irc)
  if (irc /= 0) return  ! exit (decomposition error)

  ! Check diagonal and set up error statistics.
  ! -------------------------------------------

  call CD_Diag(CD_Vec,.true.,Converged,Thr,ThrNeg,ThrFail,Diag,Qual(:,0),Buf,nDim,lBuf,ErrStat,NumCho,irc)
  if (irc /= 0) then
    irc = irc+200 ! makes the error code 400 + n
    return ! exit (check error)
  else if (.not. Converged) then
    irc = 1
    return ! exit (decomposition failure)
  end if

end if

! That's it!
! ----------

end subroutine ChoDec
