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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "numbers.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
dimension SmatPure(*)
logical InCutOff

! Enter.

! Common section.

Scalar = 0

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
        kaunter = i*(i+1)/2
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
      if (i >= j) then
        kaunter = i*(i+1)/2-i+j
      else
        kaunter = j*(j+1)/2-j+i
      end if
      ind1 = nSize*(nEqState-1)+i-1
      ind2 = nSize*(nEqState-1)+j-1
      Scalar = Scalar+Work(iVecs+ind1)*Work(iVecs+ind2)*SmatPure(kaunter)
    end do
  end do
  AddRep = exrep4*abs(Scalar)**2+exrep6*abs(Scalar)**3+exrep10*abs(Scalar)**5
end if


! Crazy energy added if inner cut-off has been passed. Ensure reject.

if (InCutOff) AddRep = 1D+20

return

end subroutine BoostRep
