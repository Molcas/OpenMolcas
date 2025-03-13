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
! Copyright (C) 2024, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  TrGrad
!
!> @brief Transform a gradient to a rotated reference, all symmetry blocks
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Transforms a gradient computed with respect to the rotation parameters
!> at the "current" orbitals to the rotation parameters at the "reference"
!> orbitals.
!>
!> @param[in]    kapOV  Parameters of the antisymmetric matrix
!> @param[in]    nkapOV Number of elements in kapOV
!> @param[inout] Grad   Gradient to transform
!> @param[in]    nOcc   Number of occupied orbitals (not including frozen) in each symmetry
!***********************************************************************

subroutine TrGrad(kapOV,nKapOV,Grad,nOcc)

use InfSCF, only: nFro, nOrb, nSym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nKapOV, nOcc(8)
real(kind=wp), intent(in) :: kapOV(nKapOV)
real(kind=wp), intent(inout) :: Grad(nKapOV)
integer(kind=iwp) :: iOff, iSym, nVir

iOff = 1

do iSym=1,nSym
  nVir = nOrb(iSym)-nFro(iSym)-nOcc(iSym)

  if (nVir*nOcc(iSym) == 0) cycle

# ifdef _SVD_
  call trg_svd(nVir,nOcc(iSym),kapOV(iOff),Grad(iOff))
# else
  call trg_series(nVir,nOcc(iSym),kapOV(iOff),Grad(iOff))
# endif

  iOff = iOff+nOcc(iSym)*nVir
end do

end subroutine TrGrad
