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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_VecBuf_CompareNormAndSum(n,nVec,Vec,J1,iSym,irc)
!
! Thomas Bondo Pedersen, September 2012.
!
! Compare norm and sum of vectors against values stored in buffer
! for Cholesky vectors J1,J1+1,J1+nVec-1 of symmetry iSym stored in
! array Vec(n,nVec). Comparison is only made for those vectors that
! are actually stored in the buffer, of course.
!
! Return codes:
!    irc=0: no differences detected.
!    irc>0: number of vectors for which differences are detected.

use ChoVecBuf, only: CHVBFI, ip_CHVBFI_SYM, nVec_in_Buf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n, nVec, J1, iSym, irc
real(kind=wp) :: Vec(n,nVec)
integer(kind=iwp) :: J, J0, mVec
real(kind=wp) :: Nrm, Sm
real(kind=wp), parameter :: Tol = 1.0e-12_wp
real(kind=wp), external :: Cho_dSumElm, dDot_
! Statement functions
real(kind=wp) :: RefNorm, RefSum
integer(kind=iwp) :: k, l
RefNorm(k,l) = CHVBFI(ip_ChVBfI_Sym(l)+2*(k-1))
RefSum(k,l) = CHVBFI(ip_ChVBfI_Sym(l)+2*(k-1)+1)

irc = 0
if (allocated(CHVBFI)) then
  J0 = J1-1
  mVec = min(J0+nVec,nVec_in_Buf(iSym))-J0
  do J=1,mVec
    Nrm = sqrt(dDot_(n,Vec(1,J),1,Vec(1,J),1))
    Sm = Cho_dSumElm(Vec(1,J),n)
    if ((abs(RefNorm(J0+J,iSym)-Nrm) > Tol) .or. (abs(RefSum(J0+J,iSym)-Sm) > Tol)) irc = irc+1
  end do
end if

end subroutine Cho_VecBuf_CompareNormAndSum
