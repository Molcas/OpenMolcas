!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine BoostRep(AddRep,SmatPure,Vecs,nSize,InCutOff)

use qmstat_global, only: exrep10, exrep4, exrep6, iOcc1, nEqState, QmType
use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: AddRep
integer(kind=iwp), intent(in) :: nSize
real(kind=wp), intent(in) :: SmatPure(*), Vecs(nSize,nSize)
logical(kind=iwp), intent(in) :: InCutOff
integer(kind=iwp) :: i, iO1, iO2, j, kaunter
real(kind=wp) :: Scalar

! Enter.

! Common section.

Scalar = Zero

! Take different route for different QM-method.

if (QmType(1:3) == 'SCF') then

  ! Repulsion term added to the Energy
  ! Calculated as S_i*S_i*CPsi_m* CPis_n
  ! S_i is the overlap integral (with the solvent molecule)
  ! for the occupied orbitals of the quantum system
  ! CPsi_m and CPsi_n are the transformation coefficients
  ! obtained from the diagonalization procedure of the Fock matrix
  ! to go from the original wavefunction to the final wavefunction
  ! after the SCF procedure. These coeficientes run over all basis set

  do iO1=1,nSize
    do iO2=1,nSize
      do i=1,iOcc1
        kaunter = nTri_Elem(i)
        Scalar = Scalar+(Vecs(i,iO1)*Vecs(i,iO2)*SmatPure(kaunter))
      end do
    end do
  end do
  Scalar = abs(Scalar)
  AddRep = exrep4*Scalar**2+exrep6*Scalar**3+exrep10*Scalar**5
else if (QmType(1:4) == 'RASS') then
  do i=1,nSize
    do j=1,nSize
      kaunter = iTri(i,j)
      Scalar = Scalar+Vecs(i,nEqState)*Vecs(j,nEqState)*SmatPure(kaunter)
    end do
  end do
  Scalar = abs(Scalar)
  AddRep = exrep4*Scalar**2+exrep6*Scalar**3+exrep10*Scalar**5
end if

! Crazy energy added if inner cut-off has been passed. Ensure reject.

if (InCutOff) AddRep = huge(AddRep)

return

end subroutine BoostRep
