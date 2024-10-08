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
!> @param[in]  mynOcc Number of occupied orbitals (including frozen) in each symmetry
!***********************************************************************

#define qnext
subroutine ExpKap(kapOV,nKapOV,U,mynOcc)

use InfSCF, only: nFro, nOFs, nOrb, nSym, TimFld
use Constants, only: Zero, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nKapOV, mynOcc(8)
real(kind=wp), intent(in) :: kapOV(nkapOV)
real(kind=wp), intent(out) :: U(nOFS)
integer(kind=iwp) :: iKap, iSym, iU, j, jU, mOcc, mOrb, mVir
real(kind=wp) :: Cpu1, Cpu2, Tim1, Tim2, Tim3
#ifndef qnext
real(kind=wp) :: theta
#endif
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

do j=1,nKapOV
  if (abs(KapOV(j)) > Pi) then
    write(u6,*) 'ExpKap: KapOV too large:',KapOV(j)
    call Abend()
  end if
end do
call Timing(Cpu1,Tim1,Tim2,Tim3)

iU = 1
iKap = 1
U(:) = Zero

do iSym=1,nSym
  mOrb = nOrb(iSym)-nFro(iSym)
  mOcc = mynOcc(iSym)-nFro(iSym)
  mVir = mOrb-mOcc

  if (mVir*mOcc == 0) cycle

  jU = iU+mOcc

  do j=1,mOcc
    U(jU:jU+mVir-1) = kapOV(iKap:iKap+mVir-1)
    iKap = iKap+mVir
    jU = jU+mOrb
  end do

# ifdef  qnext
  call matexp(mOrb,mOcc,U(iU:iU+mOrb**2))
# else
  call Exp_Schur(mOrb,U(iU:iU+mOrb**2),theta)
# endif

  iU = iU+mOrb**2
end do

do j=1,nOFS
  if (abs(U(j)) < Thrs) U(j) = Zero
end do

call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(10) = TimFld(10)+(Cpu2-Cpu1)

return

end subroutine ExpKap
