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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
!  CD_InCore_p
!
!> @brief
!>   Cholesky decompose the symmetric positive (semi-)definite matrix \p X
!> @author Francesco Aquilante
!>
!> @details
!> The \p n &times; \p n matrix \p X is Cholesky decomposed and the resulting
!> \p NumCho Cholesky vectors are returned in array \p Vec. Note that \p X
!> is modified in this routine! The array \p iD(k) on exit contains
!> the index of the diagonal from which the \p k -th Cholesky
!> vector was generated. A non-zero return code signals
!> that an error has occured (\p X might e.g. be non-positive
!> (semi-)definite) and the output is ill-defined.
!>
!> @param[in,out] X      Matrix to be Cholesky decomposed
!> @param[in]     n      Linear dimension of \p X
!> @param[out]    Vec    Storage array for Cholesky vectors
!> @param[in]     MxVec  Max. number of vectors allowed
!> @param[in,out] iD     Index array for parent diagonals
!> @param[out]    NumCho Number of Cholesky vectors
!> @param[in]     Thr    Decomposition threshold (precision)
!> @param[out]    irc    Return code
!***********************************************************************

subroutine CD_InCore_p(X,n,Vec,MxVec,iD,NumCho,Thr,irc)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, MxVec
real(kind=wp), intent(inout) :: X(n,n)
real(kind=wp), intent(out) :: Vec(n,MxVec)
integer(kind=iwp), intent(inout) :: iD(MxVec)
integer(kind=iwp), intent(out) :: NumCho, irc
real(kind=wp), intent(in) :: Thr
real(kind=wp) :: Thr_
real(kind=wp), parameter :: DefThr = 1.0e-6_wp, ThrFail = -1.0e-8_wp, ThrNeg = -1.0e-13_wp

irc = 0
NumCho = 0
if (n >= 1) then
  Thr_ = Thr
  if (Thr_ < Zero) Thr_ = DefThr

  if (MxVec > 0) then
    call CD_InCore_1p(X,n,Vec,MxVec,NumCho,Thr_,ThrNeg,ThrFail,iD,irc)
  else
    irc = -1
  end if
end if

end subroutine CD_InCore_p
