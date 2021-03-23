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

implicit real*8(a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
real*8 Dmat(*)
integer nBas
character*(LENIN4) Name(mxBas)
integer, allocatable :: nBas_per_Atom(:), nBas_Start(:)
real*8, allocatable :: Charge_B(:)
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_iScalar('Unique atoms',nAtoms)

if (nAtoms < 1) then
  write(6,'(A,I9)') 'nUniqAt =',nAtoms
  call Abend()
end if

call mma_allocate(nBas_per_Atom,nAtoms,Label='nBpA')
call mma_allocate(nBas_Start,nAtoms,Label='nB_Start')

call Get_cArray('Unique Basis Names',Name,(LENIN8)*nBas)

call BasFun_Atom(nBas_per_Atom,nBas_Start,Name,nBas,nAtoms,.false.)

call mma_allocate(Charge_B,nAtoms,Label='Charge_B')
call Get_dArray('Nuclear charge',Charge_B,nAtoms)

ZA = 0.0d0
iAt = 1
do while ((iAt <= nAtoms) .and. (ZA == 0.0d0))
  ZA = Charge_B(iAt)
  iAt = iAt+1
end do
iAt_B = iAt-1  ! start of atoms of subsystem B
call mma_deallocate(Charge_B)

if (iAt_B == 1) then ! subsystem B comes first
  ZB = 1.0d0
  nAt_B = 1
  do while ((nAt_B <= nAtoms) .and. (ZB > 0.0d0))
    ZB = Charge_B(nAt_B)
    nAt_B = nAt_B+1
  end do
  nAt_B = nAt_B-1  ! end of atoms of subsystem B
  nBas_B = nBas_Start(nAt_B)-1
  do j=nBas_B,nBas-1
    jj = j*(j+1)/2
    do i=1,j
      ijj = i+jj
      Dmat(ijj) = 0.0d0
    end do
  end do
else
  nBas_A = nBas_Start(iAt_B)-1
  nAA = nBas_A*(nBas_A+1)/2
  call FZero(Dmat,nAA)
  do j=nBas_A,nBas-1
    jj = j*(j+1)/2
    do i=1,nBas_A
      ijj = i+jj
      Dmat(ijj) = 0.0d0
    end do
  end do
end if

call mma_deallocate(nBas_Start)
call mma_deallocate(nBas_per_Atom)

! Annihilated density written to runfile for use in Coulomb gradients

Length = nBas*(nBas+1)/2
call Put_D1ao_Var(Dmat,Length)

return

end subroutine Annihil_rho
