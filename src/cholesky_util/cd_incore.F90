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
!  CD_InCore
!
!> @brief
!>   Cholesky decompose the symmetric positive (semi-)definite matrix \p X
!> @author Thomas Bondo Pedersen
!>
!> @details
!> The \p n &times; \p n matrix \p X is Cholesky decomposed and the resulting
!> \p NumCho Cholesky vectors are returned in array \p Vec. Note that \p X
!> is modified in this routine! A non-zero return code signals
!> that an error has occured (\p X might e.g. be non-positive
!> (semi-)definite) and the output is ill-defined.
!>
!> @param[in,out] X      Matrix to be Cholesky decomposed
!> @param[in]     n      Linear dimension of \p X
!> @param[out]    Vec    Storage array for Cholesky vectors
!> @param[in]     MxVec  Max. number of vectors allowed
!> @param[out]    NumCho Number of Cholesky vectors
!> @param[in]     Thr    Decomposition threshold (precision)
!> @param[out]    irc    Return code
!***********************************************************************

subroutine CD_InCore(X,n,Vec,MxVec,NumCho,Thr,irc)

implicit none
integer n, MxVec, NumCho, irc
real*8 X(n,n)
real*8 Vec(n,MxVec)
real*8 Thr
character*9 SecNam
parameter(SecNam='CD_InCore')
real*8 DefThr
parameter(DefThr=1.0d-6)
real*8 ThrNeg, ThrFail
parameter(ThrNeg=-1.0d-13,ThrFail=-1.0d-8)

irc = 0
NumCho = 0
if (n < 1) Go To 1 ! exit (nothing to do)
if (Thr < 0.0d0) Thr = DefThr

if (MxVec > 0) then
  call CD_InCore_1(X,n,Vec,MxVec,NumCho,Thr,ThrNeg,ThrFail,irc)
else
  irc = -1
end if

1 continue

end subroutine CD_InCore
