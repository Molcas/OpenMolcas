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
subroutine rhodyn()
!***********************************************************************
!     module rhodyn_data contains common variables declaration
!     and also comments on their meaning
!***********************************************************************
  use rhodyn_data
  use rhodyn_utils, only: mult, dashes, sortci
  use definitions, only: wp, iwp, u6
  use constants, only: auToeV
  use stdalloc, only: mma_allocate, mma_deallocate
  use mh5, only: mh5_put_dset, mh5_close_file
  implicit none
  integer(kind=iwp) :: lu !temporary io unit
  integer(kind=iwp), external :: iPrintLevel, isFreeUnit
  real(kind=wp), dimension(:,:), allocatable :: Ham ! auxiliary matrix

!  call qEnter('RHODYN') !deprecated
  ireturn = 20
  ipglob  = iPrintLevel(-1)

  call StatusLine('RHODYN:','Starting calculation')

! initializitation of default values and printing main parameters
  call rhodyn_init()
! reading input file
  call read_input()

! calculation of integer parameters needed for allocation
  maxnum   = maxval(ndet)
  maxnconf = maxval(nconf)
  maxlroots= maxval(lroots)
  nconftot =0
  lrootstot=0
  if (p_style=='DET') ndet_tot=sum(NDET)   ! Nr of DETs
  do i=1,N
    if (flag_so) then
      nconftot=nconftot+nconf(i)*ispin(i)    ! Nr of CSFs
      lrootstot=lrootstot+lroots(i)*ispin(i) ! Nr of SF States
    else
      lrootstot=lrootstot+lroots(i)
      nconftot=nconftot+nconf(i)
    endif
  enddo

  call mma_allocate(dipole, lrootstot, lrootstot,3)
  call mma_allocate(dysamp, lrootstot, lrootstot)
  call mma_allocate(dysamp_bas, lrootstot, lrootstot)
  call mma_allocate(tmp,       Nstate,   Nstate     )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     start from rassf/rassi output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (preparation/=2.and.preparation/=4) then

    n_freq = Nstate*(Nstate-1)/2

    call mma_allocate(H_CSF, maxnconf, maxnconf, N)
    call mma_allocate(CI, maxnconf, maxlroots, N)
    call mma_allocate(E, maxlroots,N)
    call mma_allocate(U_CI, nconftot, lrootstot)
    call mma_allocate(SO_CI, lrootstot, lrootstot)
    call mma_allocate(DTOC, maxnum, maxnconf, N)
    call mma_allocate(V_CSF, nconftot, nconftot)
    call mma_allocate(V_SO, lrootstot, lrootstot)
    call mma_allocate(CSF2SO, nconftot, lrootstot)
    call mma_allocate(E_SO, lrootstot)
    call mma_allocate(HTOTRE_CSF, nconftot, nconftot)
    call mma_allocate(HTOT_CSF, nconftot, nconftot)
    call mma_allocate(a_einstein, lrootstot, lrootstot)
    call mma_allocate(emiss, n_freq)

    call cre_prep()

    H_CSF=0d0
    CI=0d0

! Determine file names
! Expected N rassd files and 1 rassisd file
    call mma_allocate(rassd_list,N)
    do i=1,N
      write(rassd_list(i),"(A5,I1)") "RASSD", i
      write(u6,*) rassd_list(i)
! c
      write(u6,*) 'entering read_rassd: ', i
      call read_rassd(i)
    enddo
    call mma_deallocate(rassd_list)

    if (i_rasscf==1) then
! currently disabled
! Diagonalize H(RASSCF) and sort ascending the eigenvalue get the
! corresponding eigenvector (CI coefficients) to matrix CI(:,:,N)
      do i=1,N
        if (ipglob>3) then
          write(u6,sint)'H(RASSCF) in CSF basis of spin manifold:',i
          call dashes()
          do k=1,15
            write(u6,*)(H_CSF(k,j,i),j=1,15)
          enddo
          call dashes()
        endif
        call mma_allocate(Ham,lroots(i),lroots(i))
        Ham(:,:)=H_CSF(1:nconf(i),1:nconf(i),i)
        call sortci(nconf(i), Ham, E(1:nconf(i),i), &
                            CI(1:nconf(i),1:lroots(i),i),ipglob)
        call mma_deallocate(Ham)
      enddo
    else if (i_rasscf==2.or.i_rasscf==3) then
! Construct SF Hamiltonians from CIs and Es for each spin manifold
      do i=1,N
        do j=1,lroots(i)
          H_CSF(j,j,i) = E(j,i)
        enddo
! H_CSF=CI*diag(E)*CI^T
        call mma_allocate(Ham,nconf(i),lroots(i))
        call mult(CI(1:nconf(i),1:lroots(i),i), &
                  H_CSF(1:lroots(i),1:lroots(i),i),Ham)
        call mult(Ham,CI(1:nconf(i),1:lroots(i),i), &
                  H_CSF(1:nconf(i),1:nconf(i),i),.False.,.True.)
        call mma_deallocate(Ham)
      enddo

    endif

    call mh5_put_dset(prep_ci, CI)
    call mh5_put_dset(prep_hcsf, H_CSF)

    call uci()
    if (flag_so) call transform_V()
    call get_hcsf()
    if (flag_so) then
      call soci()
      call get_dipole()
    endif
    call get_dm0()

    call mh5_close_file(prep_id)

    call dashes()
    call dashes()
    write(u6,*) 'Preparation finished successfully'
    call dashes()
    call dashes()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start from intermadiate preparation file PREP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else if (preparation==2) then
    call mma_allocate(HTOT_CSF, nconftot, nconftot)
    call read_prep()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! charge migration case
! only SO hamiltonian from RASSI is read
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else if (preparation==4) then
! only SO basis allowed
  !basis = 'SO'
    call mma_allocate(HSOCX, Nstate, Nstate)
    call mma_allocate(E_SO, Nstate)
    call get_dipole()
    if (DM_basis=='CSF_SO') then
      call mma_allocate(CI, maxnconf, maxlroots,N)
      call mma_allocate(E, maxlroots, N)
      call mma_allocate(U_CI, nconftot, lrootstot)
      call mma_allocate(CSF2SO, nconftot, lrootstot)
      CI=0d0
      ! Determine file names
      ! Expected N rassd files and 1 rassisd file
      call mma_allocate(rassd_list,N)
      do i=1,N
        write(rassd_list(i),"(A5,I1)") "RASSD", i
        write(u6,*) rassd_list(i)
        ! c
        write(u6,*) 'entering read_rassd: ', i
        call read_rassd(i)
      enddo
      call mma_deallocate(rassd_list)
      ! construct CI transformation matrix
      call uci()
      call read_rassisd()
      ! construct the transformation matrix CSF2SO from SO to CSF
      call mult(dcmplx(U_CI),SO_CI,CSF2SO)
    endif
!  call read_rassisd()
    call get_dm0()
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dynamics part starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (preparation/=3) then
    select case (basis)
    case ('CSF')
      d = nconftot
    case ('SF')
      d = lrootstot
    case ('SO')
      d = lrootstot
    end select
    call mma_allocate(hamiltonian, d,d)
    call mma_allocate(density0,    d,d)
    call mma_allocate(dipole_basis,d,d,3)
    call mma_allocate(hamiltoniant,Nstate,Nstate)
    call mma_allocate(densityt,    Nstate,Nstate)
    call mma_allocate(decay,       Nstate,Nstate)
    call mma_allocate(U_CI_compl,  nconftot,lrootstot)

! prepare density and hamiltonian in required basis to propagate with
    if (preparation/=4) then
      call hamdens()
    else
    ! charge migration case
      hamiltonian = HSOCX
      density0 = DM0
      dipole_basis= dipole
      U_CI_compl =dcmplx(U_CI,0d0)
      if (flag_dyson) then
        do i=1,3
          dipole_basis(:,:,i)=dipole_basis(:,:,i)+alpha*dysamp_bas
        enddo
      endif
    endif
! dynamics will be performed internally using matrices
! of dimension d
    if (basis=='CSF') then
      ! call assert(nconftot==Nstate)
      d = nconftot
    elseif (basis=='SF'.or.basis=='SO') then
      if (Nstate==lrootstot) then
        d = lrootstot
      elseif (Nstate<lrootstot) then
        d = Nstate
        call cut_matrices()
      endif
    endif

    if (flag_diss) then
      call mma_allocate(kab_basis,d,d)
      call mma_allocate(k_bar_basis,d,d)
      if (kext) then
        call k_external()
      else
        call kab()
      endif
      if (ipglob>3) then
        lu = isFreeUnit(19)
        call molcas_open(lu,'Kab_basis_eV.dat')
        do i=1,d
          write(lu,*)(kab_basis(i,j)*autoev,j=1,d)
        enddo
        close(lu)
        lu = isFreeUnit(20)
        call molcas_open(lu,'K_bar_basis_eV.dat')
        do i=1,d
          write(lu,*)(k_bar_basis(i,j)*autoev,j=1,d)
        enddo
        close(lu)
      endif
    endif

    if (flag_decay.or.flag_dyson) call prepare_decay()

    call propagate()

  endif

! closing and deallocation

  call mh5_close_file(out_id)

! allocated in read_input
  if (allocated(ndet)) call mma_deallocate(ndet)
  if (allocated(nconf)) call mma_deallocate(nconf)
  if (allocated(lroots)) call mma_deallocate(lroots)
  if (allocated(ispin)) call mma_deallocate(ispin)
  if (allocated(istates)) call mma_deallocate(istates)
  if (allocated(taushift)) call mma_deallocate(taushift)
  if (allocated(amp)) call mma_deallocate(amp)
  if (allocated(sigma)) call mma_deallocate(sigma)
  if (allocated(omega)) call mma_deallocate(omega)
  if (allocated(phi)) call mma_deallocate(phi)
  if (allocated(pulse_vector)) call mma_deallocate(pulse_vector)
!
! allocated in prepare part
  if (allocated(V_CSF)) call mma_deallocate(V_CSF)
  if (allocated(V_SO)) call mma_deallocate(V_SO)
  if (allocated(H_CSF)) call mma_deallocate(H_CSF)
  if (allocated(CI)) call mma_deallocate(CI)
  if (allocated(DTOC)) call mma_deallocate(DTOC)
  if (allocated(E)) call mma_deallocate(E)
  if (allocated(U_CI)) call mma_deallocate(U_CI)
  if (allocated(SO_CI)) call mma_deallocate(SO_CI)
  if (allocated(E_SO)) call mma_deallocate(E_SO)
  if (allocated(HTOTRE_CSF)) call mma_deallocate(HTOTRE_CSF)
  if (allocated(dipole)) call mma_deallocate(dipole)
  if (allocated(dysamp)) call mma_deallocate(dysamp)
  if (allocated(dysamp_bas)) call mma_deallocate(dysamp_bas)
  if (allocated(HTOT_CSF)) call mma_deallocate(HTOT_CSF)
  if (allocated(CSF2SO)) call mma_deallocate(CSF2SO)
  if (allocated(tmp)) call mma_deallocate(tmp)
  if (allocated(DM0)) call mma_deallocate(DM0)
  if (allocated(HSOCX)) call mma_deallocate(HSOCX)

! allocated in dynamics part
  if (preparation/=3) then
    if (allocated(U_CI_compl)) call mma_deallocate(U_CI_compl)
    if (allocated(dipole)) call mma_deallocate(dipole)
    if (allocated(hamiltonian)) call mma_deallocate(hamiltonian)
    if (allocated(hamiltoniant)) call mma_deallocate(hamiltoniant)
    if (allocated(density0)) call mma_deallocate(density0)
    if (allocated(decay)) call mma_deallocate(decay)
    if (allocated(densityt)) call mma_deallocate(densityt)
    if (allocated(dipole_basis)) call mma_deallocate(dipole_basis)
    if (allocated(a_einstein)) call mma_deallocate(a_einstein)
    if (allocated(emiss)) call mma_deallocate(emiss)
    if (allocated(kab_basis)) call mma_deallocate(kab_basis)
    if (allocated(k_bar_basis)) call mma_deallocate(k_bar_basis)
  endif


  call StatusLine('RhoDyn:','Finished')
  ireturn = 0
!  call qExit('RHODYN') ! deprecated
  return
end
