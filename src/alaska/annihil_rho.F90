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

subroutine Annihil_rho(Dmat,nBas)

#include "intent.fh"

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(_OUT_) :: Dmat(*)
integer(kind=iwp), intent(in) :: nBas
#include "Molcas.fh"
integer(kind=iwp) :: i, iAt, iAt_B, ijj, j, jj, Length, nAA, nAt_B, nAtoms, nBas_A, nBas_B
real(kind=wp) :: ZA, ZB
integer(kind=iwp), allocatable :: nBas_per_Atom(:), nBas_Start(:)
real(kind=wp), allocatable :: Charge_B(:)
character(len=LenIn8), allocatable :: UBName(:)

!                                                                      *
!***********************************************************************
!                                                                      *
call Get_iScalar('Unique atoms',nAtoms)

if (nAtoms < 1) then
  write(u6,'(A,I9)') 'nUniqAt =',nAtoms
  call Abend()
end if

call mma_allocate(nBas_per_Atom,nAtoms,Label='nBpA')
call mma_allocate(nBas_Start,nAtoms,Label='nB_Start')

call mma_allocate(UBName,nBas,Label='UBName')
call Get_cArray('Unique Basis Names',UBName,(LENIN8)*nBas)

call BasFun_Atom(nBas_per_Atom,nBas_Start,UBName,nBas,nAtoms,.false.)
call mma_deallocate(UBName)

call mma_allocate(Charge_B,nAtoms,Label='Charge_B')
call Get_dArray('Nuclear charge',Charge_B,nAtoms)

ZA = Zero
iAt = 1
do while ((iAt <= nAtoms) .and. (ZA == Zero))
  ZA = Charge_B(iAt)
  iAt = iAt+1
end do
iAt_B = iAt-1  ! start of atoms of subsystem B
call mma_deallocate(Charge_B)

if (iAt_B == 1) then ! subsystem B comes first
  ZB = One
  nAt_B = 1
  do while ((nAt_B <= nAtoms) .and. (ZB > Zero))
    ZB = Charge_B(nAt_B)
    nAt_B = nAt_B+1
  end do
  nAt_B = nAt_B-1  ! end of atoms of subsystem B
  nBas_B = nBas_Start(nAt_B)-1
  do j=nBas_B,nBas-1
    jj = j*(j+1)/2
    do i=1,j
      ijj = i+jj
      Dmat(ijj) = Zero
    end do
  end do
else
  nBas_A = nBas_Start(iAt_B)-1
  nAA = nBas_A*(nBas_A+1)/2
  DMat(1:nAA) = Zero
  do j=nBas_A,nBas-1
    jj = j*(j+1)/2
    do i=1,nBas_A
      ijj = i+jj
      Dmat(ijj) = Zero
    end do
  end do
end if

call mma_deallocate(nBas_Start)
call mma_deallocate(nBas_per_Atom)

! Annihilated density written to runfile for use in Coulomb gradients

Length = nBas*(nBas+1)/2
call Put_dArray('D1aoVar',Dmat,Length)

return

end subroutine Annihil_rho
