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
! Copyright (C) Ben Swerts                                             *
!               2016, Liviu Ungur                                      *
!***********************************************************************

subroutine MakeDens(nBas,nOrb,Cff,OrbEn,EnergyWeight,Dens)
!***********************************************************************
!                                                                      *
!     purpose: Compute (energy weighted) density matrix in AO basis    *
!                                                                      *
!     input:                                                           *
!       nBas    : number of basis functions                            *
!       nOrb    : number of orbitals                                   *
!       Cff     : molecular orbital coefficients                       *
!       OrbEn   : molecular orbital energies                           *
!       EnergyWeight : if .true. do energy weighting                   *
!                                                                      *
!     output:                                                          *
!       Dens    : density matrix in triangular storage                 *
!                                                                      *
!     called from: Drv2El_FAIEMP                                       *
!                  FragPInt                                            *
!                                                                      *
!     written by: B. Swerts                                            *
!     modified by L. Ungur                                             *
!     simplified version of scf/done_scf.f                             *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas, nOrb
real(kind=wp), intent(in) :: Cff(nBas*nOrb), OrbEn(nOrb)
logical(kind=iwp), intent(in) :: EnergyWeight
real(kind=wp), intent(out) :: Dens(nBas*(nBas+1)/2)
integer(kind=iwp) :: i, Ind, ij, iRow, iCol
real(kind=wp) :: energy, rSum
logical(kind=iwp) :: DBG

energy = One

DBG = .false.
if (DBG) write(u6,*) 'MakeDens:  EnergyWeight',EnergyWeight
if (DBG) call xFlush(u6)
if (DBG) write(u6,*) 'MakeDens:  nBas=',nBas
if (DBG) call xFlush(u6)
if (DBG) write(u6,*) 'MakeDens:  nOrb=',nOrb
if (DBG) call xFlush(u6)
if (DBG) write(u6,*) 'MakeDens: OrbEn=',(OrbEn(i),i=1,nOrb)
if (DBG) call xFlush(u6)
if (DBG) write(u6,*) 'MakeDens:   Cff=',(Cff(i),i=1,nBas*nOrb)
if (DBG) call xFlush(u6)

do iRow=1,nBas
  rSum = Zero
  ij = -1
  do i=1,nOrb
    ij = ij+1
    if (EnergyWeight) energy = OrbEn(i)
    rSum = rSum+energy*Cff(iRow+ij*nBas)*Cff(iRow+ij*nBas)
  end do
  Ind = iRow*(iRow-1)/2
  Dens(Ind+iRow) = Two*rSum

  do iCol=1,iRow-1
    rSum = Zero
    ij = -1
    do i=1,nOrb
      ij = ij+1
      if (EnergyWeight) energy = OrbEn(i)
      rSum = rSum+energy*Cff(iRow+ij*nBas)*Cff(iCol+ij*nBas)
    end do
    Dens(Ind+iCol) = Four*rSum
  end do
end do
if (DBG) call TriPrt('Dens in MakeDens',' ',Dens,nBas)
if (DBG) call xFlush(u6)

return

end subroutine MakeDens
