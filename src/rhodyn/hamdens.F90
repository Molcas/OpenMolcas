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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************

subroutine hamdens()
! Purpose : Transform initial Hamiltonian to required basis as well as
! the initial density matrix, which are currently present
! in CSF basis

use rhodyn_data, only: alpha, basis, CSF2SO, d, density0, dipole, dipole_basis, DM0, dysamp_bas, flag_dyson, flag_so, &
                       hamiltonian, HTOT_CSF, initialtime, ipglob, Nstate, SO_CI, tmp, U_CI, U_CI_compl
use rhodyn_utils, only: transform, dashes
use linalg_mod, only: mult
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, j, k, ii

U_CI_compl(:,:) = cmplx(U_CI,kind=wp)
hamiltonian(:,:) = cZero
density0(:,:) = cZero

! construct the initial hamiltonian and density matrix

write(u6,*) 'Basis: ',basis
if (initialtime == Zero) then
  if (basis == 'CSF') then
    hamiltonian(:,:) = HTOT_CSF
    density0(:,:) = DM0
  else if (basis == 'SF') then
    ! Hamiltonian CSF->SF
    call transform(HTOT_CSF,U_CI_compl,hamiltonian)
    ! density CSF->SF
    call transform(DM0,U_CI_compl,density0)
  else if (basis == 'SO') then
    ! Hamiltonian CSF->SO
    call transform(HTOT_CSF,CSF2SO,hamiltonian)
    ! density CSF->SO
    call transform(DM0,CSF2SO,density0)
  end if
else if (initialtime /= Zero) then
  ! if initialtime /= 0, then the initial density matrix in basis
  ! of "Basis", not in CSF basis, in this case, the density matrix
  ! readin was printout directly from the propagate part, so we can
  ! directly use it in propagate section, let density0 equal to the read
  ! in density matrix.
  if (basis == 'CSF') then
    ! Hamiltonian CSF
    hamiltonian(:,:) = HTOT_CSF
    ! density SF->CSF
    call transform(DM0,U_CI_compl,density0,.false.)
  else if (basis == 'SF') then
    ! Hamiltonian CSF->SF
    call transform(HTOT_CSF,U_CI_compl,hamiltonian)
    ! density SF
    density0(:,:) = DM0
  else if (basis == 'SO') then
    ! Hamiltonian CSF->SO
    call transform(HTOT_CSF,CSF2SO,hamiltonian)
    ! density SO
    density0(:,:) = DM0
  end if
end if

! construct the dipole matrix in required basis (with dyson matrix)
! here dipole presumably is in SO basis (if SO is on) or in SF basis (is SO is off)

if (flag_so) then
  if (basis == 'CSF') then
    do i=1,3
      call transform(dipole(:,:,i),CSF2SO,dipole_basis(:,:,i),.false.)
    end do
  else if (basis == 'SF') then
    do i=1,3
      call transform(dipole(:,:,i),SO_CI,dipole_basis(:,:,i),.false.)
    end do
  else if (basis == 'SO') then
    dipole_basis(:,:,:) = dipole
    if (flag_dyson) then
      call mult(SO_CI,dysamp_bas,tmp,.true.,.false.)
      call mult(tmp,SO_CI,dysamp_bas)
    end if
  end if
else ! flag_so is off
  if (basis == 'CSF') then
    do i=1,3
      call transform(dipole(:,:,i),U_CI_compl,dipole_basis(:,:,i),.false.)
    end do
  else if (basis == 'SF') then
    dipole_basis(:,:,:) = dipole
  end if
end if

if (flag_dyson) then
  do i=1,3
    dipole_basis(:,:,i) = dipole_basis(:,:,i)+alpha*dysamp_bas
  end do
end if

if (ipglob > 3) then
  ii = 10
  if (Nstate < 10) ii = Nstate
  write(u6,*) 'hamiltonian'
  do i=1,ii
    write(u6,*) (hamiltonian(i,j),j=1,ii)
  end do
  write(u6,*) 'density0'
  do i=1,ii
    write(u6,*) (density0(i,j),j=1,ii)
  end do
  write(u6,*) 'End get_dipole'
  call dashes()
end if
if (ipglob > 4) then
  do i=1,3
    call dashes()
    write(u6,*) 'Dipole Matrix in',basis,'basis'
    if (i == 1) write(u6,*) 'Printout the components dipole matrix dx'
    if (i == 2) write(u6,*) 'Printout the components dipole matrix dy'
    if (i == 3) write(u6,*) 'Printout the components dipole matrix dz'
    call dashes()
    do k=1,d
      write(u6,*) (dipole_basis(k,j,i),j=1,d)
    end do
    write(u6,*)
    call dashes()
  end do
end if

end subroutine hamdens
