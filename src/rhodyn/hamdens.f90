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
  use rhodyn_data
  use rhodyn_utils, only: transform, mult, dashes
  use definitions, only: u6
  implicit none
!
! Purpose : Transform initial Hamiltonian to required basis as well as
! the initial density matrix, which are currently present
! in CSF basis
!
  if (ipglob>2) then
    call dashes()
    write(u6,*) 'Begin get_hamiltonian'
    call dashes()
  endif

  if (.not.flag_test) then
    U_CI_compl(:,:) =dcmplx(U_CI,0d0)
  endif
  hamiltonian=zero
  density0=zero

! construct the initial hamiltonian and density matrix

  write(u6,*) 'Basis: ', basis
  if (initialtime==0d0) then
    if (basis=='CSF') then
      hamiltonian(:,:) = HTOT_CSF
      density0(:,:) = DM0
    elseif (basis=='SF') then
      ! Hamiltonian CSF->SF
      call transform(HTOT_CSF,U_CI_compl,hamiltonian)
      ! density CSF->SF
      call transform(DM0,U_CI_compl,density0)
    elseif (basis=='SO') then
      ! Hamiltonian CSF->SO
      call transform(HTOT_CSF,CSF2SO,hamiltonian)
      ! density CSF->SO
      call transform(DM0,CSF2SO,density0)
    endif
  elseif (initialtime/=0d0) then
! if initialtime .ne. 0, then the initial density matrix in basis
! of "Basis", not in CSF basis, in this case, the density matrix
! readin was printout directly from the propagate part, so we can
! directly use it in propagate section, let density0 equal to the read
! in density matrix.
    if (basis=='CSF') then
      ! Hamiltonian CSF
      hamiltonian(:,:) = HTOT_CSF
      ! density SF->CSF
      call transform(DM0,U_CI_compl,density0,.False.)
    elseif (basis=='SF') then
      ! Hamiltonian CSF->SF
      call transform(HTOT_CSF,U_CI_compl,hamiltonian)
      ! density SF
      density0(:,:) = DM0
    elseif (basis=='SO') then
      ! Hamiltonian CSF->SO
      call transform(HTOT_CSF,CSF2SO,hamiltonian)
      ! density SO
      density0(:,:) = DM0
    endif
  endif

  if (ipglob>2) then
    call dashes()
    write(u6,*) 'End get_hamiltonian'
    call dashes()
  endif


  if (.not.flag_test) then

    if (ipglob>2) then
    call dashes()
    write(u6,*) 'Begin get_dipole'
    call dashes()
    endif

! contrust the dipole matrix in required basis (with dyson matrix)
! here dipole presumably in SO basis (if SO is on)
!                     or in SF basis (is SO is off)
!
    if (flag_so) then
      if (basis=='CSF') then
        do i=1,3
        call transform(dipole(:,:,i),CSF2SO,dipole_basis(:,:,i),.False.)
        enddo
      elseif (basis=='SF')then
        do i=1,3
        call transform(dipole(:,:,i),SO_CI,dipole_basis(:,:,i),.False.)
        enddo
      elseif (basis=='SO') then
        dipole_basis(:,:,:) = dipole
        if (flag_dyson) then
          call mult(SO_CI,dysamp_bas,tmp,.True.,.False.)
          call mult(tmp,SO_CI,dysamp_bas)
        endif
      endif
    else ! flag_so is off
      if (basis=='CSF') then
        do i=1,3
          call transform(dipole(:,:,i),U_CI_compl, &
                         dipole_basis(:,:,i),.False.)
        enddo
      elseif (basis=='SF') then
        dipole_basis(:,:,:) = dipole
      endif
    endif

    if (flag_dyson) then
      do i=1,3
        dipole_basis(:,:,i)=dipole_basis(:,:,i)+alpha*dysamp_bas
      enddo
    endif

!!!!!!!!!!!!!!!!!!!
    if (ipglob>3) then
      ii=10
      if (Nstate<10) ii=Nstate
      write(u6,*) 'hamiltonian'
      do i=1,ii
        write(u6,*) (hamiltonian(i,j),j=1,ii)
      enddo
      write(u6,*) 'density0'
      do i=1,ii
        write(u6,*) (density0(i,j),j=1,ii)
      enddo
      write(u6,*) 'End get_dipole'
      call dashes()
    endif
    if (ipglob>4) then
      do i=1,3
        call dashes()
        write(u6,*) 'Dipole Matrix in', basis, 'basis'
        if(i==1)write(u6,*)'Printout the components dipole matrix dx'
        if(i==2)write(u6,*)'Printout the components dipole matrix dy'
        if(i==3)write(u6,*)'Printout the components dipole matrix dz'
        call dashes()
        do k=1,d
          write(u6,*) (dipole_basis(k,j,i),j=1,d)
        enddo
        write(u6,*)
        call dashes()
      enddo
    endif
!
  endif

  if (flag_test .and. flag_pulse) then
    dipole_basis(:,:,:) = dipole
    if (ipglob>4) then
      do i=1,3
        call dashes()
        write(u6,*) 'Dipole Matrix in', basis, 'basis'
        if(i==1)write(u6,*)'Printout the components dipole matrix dx'
        if(i==2)write(u6,*)'Printout the components dipole matrix dy'
        if(i==3)write(u6,*)'Printout the components dipole matrix dz'
        call dashes()
        do k=1,d
          write(u6,*) (dipole_basis(k,j,i),j=1,d)
        enddo
        write(u6,*)
        call dashes()
      enddo
    endif
  endif

end
