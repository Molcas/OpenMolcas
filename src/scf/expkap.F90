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
! Copyright (C) 2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  ExpKap
!
!> @brief Compute an orbital rotation matrix from the rotation parameters
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Computes the orbital rotation matrix corresponding to the parametrized
!> form as an antisymmetric matrix, i.e. its exponential.
!> The input parameters \p kapOV are the unique elements of the occupied-virtual
!> block of each symmetry.
!>
!> @param[in]  kapOV  Parameters of the antisymmetric matrix
!> @param[in]  nkapOV number of elements in kapOV
!> @param[out] U      Unitary matrix to transform old CMOs
!> @param[in]  nOcc   Number of occupied orbitals (not including frozen) in each symmetry
!***********************************************************************

#define EXP_QNEXT 1
#define EXP_SERIES 2
#define EXP_SVD 3
#define EXP_FULL 4
#define _EXP_ EXP_SVD
subroutine ExpKap(kapOV,nKapOV,U,nOcc)

use InfSCF, only: nFro, nOFs, nOrb, nSym, TimFld
use Constants, only: Zero, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nKapOV, nOcc(8)
real(kind=wp), intent(in) :: kapOV(nkapOV)
real(kind=wp), intent(out) :: U(nOFS)
integer(kind=iwp) :: iKap, iSym, iU, j, jU, mOrb, mVir
logical(kind=iwp) :: use_svd
real(kind=wp) :: Cpu1, Cpu2, theta, Tim1, Tim2, Tim3
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

theta = Pi
#if ( _EXP_ == EXP_QNEXT || _EXP_ == EXP_SERIES )
! If using an expansion, fall back to SVD for large displacements
use_svd = .false.
do j=1,nKapOV
  if (abs(KapOV(j)) > Pi) then
    use_svd = .true.
    exit
  end if
end do
#elif ( _EXP_ == EXP_FULL )
use_svd = .false.
#elif ( _EXP_ == EXP_SVD )
use_svd = .true.
#endif
call Timing(Cpu1,Tim1,Tim2,Tim3)

iU = 1
iKap = 1
U(:) = Zero

do iSym=1,nSym
  mOrb = nOrb(iSym)-nFro(iSym)
  mVir = mOrb-nOcc(iSym)

  if (mVir*nOcc(iSym) == 0) cycle

  jU = iU+nOcc(iSym)

  do j=1,nOcc(iSym)
    U(jU:jU+mVir-1) = kapOV(iKap:iKap+mVir-1)
    iKap = iKap+mVir
    jU = jU+mOrb
  end do

  if (use_svd) then
    call Exp_SVD(mOrb,nOcc(iSym),U(iU:iU+mOrb**2-1),theta)
  else
#   if ( _EXP_ == EXP_QNEXT )
    call Exp_series(mOrb,nOcc(iSym),U(iU:iU+mOrb**2-1))
#   elif ( _EXP_ == EXP_SERIES )
    call Exp_series2(mOrb,nOcc(iSym),U(iU:iU+mOrb**2-1))
#   elif ( _EXP_ == EXP_FULL )
    call Exp_eig(mOrb,U(iU:iU+mOrb**2-1),theta)
#   endif
  end if

  iU = iU+mOrb**2
end do

do j=1,nOFS
  if (abs(U(j)) < Thrs) U(j) = Zero
end do

call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(10) = TimFld(10)+(Cpu2-Cpu1)

return

end subroutine ExpKap
