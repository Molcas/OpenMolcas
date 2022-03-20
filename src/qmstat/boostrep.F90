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

subroutine BoostRep(AddRep,SmatPure,iVecs,nSize,InCutOff)

use qmstat_global, only: exrep10, exrep4, exrep6, iOcc1, nEqState, QmType
use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: AddRep, SmatPure(*)
integer(kind=iwp) :: iVecs, nSize
logical(kind=iwp) :: InCutOff
#include "maxi.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ind1, ind2, iO1, iO2, j, kaunter
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
        ind1 = nSize*(iO1-1)+i-1
        ind2 = nSize*(iO2-1)+i-1
        Scalar = Scalar+(Work(iVecs+ind1)*Work(iVecs+ind2)*SmatPure(kaunter))
      end do
    end do
  end do
  AddRep = exrep4*abs(Scalar)**2+exrep6*abs(Scalar)**3+exrep10*abs(Scalar)**5
else if (QmType(1:4) == 'RASS') then
  do i=1,nSize
    do j=1,nSize
      kaunter = iTri(i,j)
      ind1 = nSize*(nEqState-1)+i-1
      ind2 = nSize*(nEqState-1)+j-1
      Scalar = Scalar+Work(iVecs+ind1)*Work(iVecs+ind2)*SmatPure(kaunter)
    end do
  end do
  AddRep = exrep4*abs(Scalar)**2+exrep6*abs(Scalar)**3+exrep10*abs(Scalar)**5
end if

! Crazy energy added if inner cut-off has been passed. Ensure reject.

if (InCutOff) AddRep = huge(AddRep)

return

end subroutine BoostRep
