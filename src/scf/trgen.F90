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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine TrGen(TrMat,nTrMat,Ovrlp,OneHam,mBT)
!***********************************************************************
!                                                                      *
!     purpose: Generate transformation matrix from AO's in which       *
!              integrals are computed to orthogonal, symmetry adapted  *
!              (and spherical, if desired) functions - near linear     *
!              dependencies are also removed.                          *
!                                                                      *
!     input:                                                           *
!       Ovrlp    : overlap matrix in AO basis of length mBT            *
!                                                                      *
!     output:                                                          *
!       TrMat   : Transformation matrix of length nTrMat               *
!                                                                      *
!***********************************************************************

use InfSCF, only: DelThr, nBas, nBO, nBT, nnFr, nSym
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTrMat, mBT
real(kind=wp), intent(out) :: TrMat(nTrMat)
real(kind=wp), intent(in) :: Ovrlp(mBT), OneHam(mBT)
integer(kind=iwp) :: i, ind, iSym, j

ind = 0
do iSym=1,nSym
  do i=1,nBas(iSym)
    do j=1,nBas(iSym)
      ind = ind+1
      TrMat(ind) = Zero
      if (I == J) TrMat(ind) = One
    end do
  end do
end do

! Set up certain parameters (nOrb(i) may be changed)
call SetUp_SCF()

! Move frozen atomic orbitals to the begining
if (nnFr > 0) then
  call Freeze(TrMat,nBO,OneHam,mBT)
  call SetUp_SCF()
end if

! Remove near linear dependencies from basis set
if (DelThr /= Zero) then
  call OvlDel(Ovrlp,nBT,TrMat,nBO)
  call SetUp_SCF()
end if

! Orthogonalize final orbitals
call Ortho(TrMat,nBO,Ovrlp,nBT)

end subroutine TrGen
