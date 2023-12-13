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
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine rhodyn(ireturn)
!***********************************************************************
! performs the calls to other subroutines that collect all data needed
! for dynamics and the second part is responsible for dynamics itself
!***********************************************************************

use rhodyn_data, only: a_einstein, alpha, amp, basis, CI, CSF2SO, d, decay, density0, densityt, dipole, dipole_basis, DM0, &
                       DM_basis, DTOC, dysamp, dysamp_bas, E, E_SF, E_SO, emiss, flag_decay, flag_diss, flag_dyson, flag_so, &
                       H_CSF, hamiltonian, hamiltoniant, HSOCX, HTOT_CSF, HTOTRE_CSF, i_rasscf, ipglob, ispin, istates, n_sf, &
                       k_bar_basis, kab_basis, kext, list_sf_mult, list_sf_spin, list_sf_states, list_so_mult, list_so_proj, &
                       list_so_sf, list_so_spin, lroots, lrootstot, maxlroots, maxnconf, N, n_freq, nconf, nconftot, ndet, &
                       ndet_tot, Nstate, omega, out_id, p_style, phi, prep_ci, prep_hcsf, prep_id, pulse_vector, rassd_list, &
                       runmode, sigma, sint, SO_CI, taushift, tmp, U_CI, U_CI_compl, V_CSF, V_SO
use rhodyn_utils, only: dashes, sortci, transform
use linalg_mod, only: mult
use mh5, only: mh5_close_file, mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, auToeV
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, ii, j, jj, k, lu, m, maxnum
! Ham: auxiliary matrix
! pop_sf: store summed populations over spin manifolds
real(kind=wp), allocatable :: Ham(:,:), pop_sf(:)
integer(kind=iwp), external :: iPrintLevel, isFreeUnit

ireturn = 20
ipglob = iPrintLevel(-1) ! MOLCAS_PRINT variable

call StatusLine('RhoDyn:','Starting calculation')

! initializitation of default values and printing main parameters
call rhodyn_init()
! reading input file
call read_input()

! calculation of parameters needed for allocation
maxnum = maxval(ndet)
maxnconf = maxval(nconf)
maxlroots = maxval(lroots)
n_sf = sum(lroots)
nconftot = 0
lrootstot = 0
if (p_style == 'DET') ndet_tot = sum(ndet)  ! Nr of DETs
do i=1,N
  if (flag_so) then
    ! total number of CSF confs and SF roots including SO splitting
    nconftot = nconftot+nconf(i)*ispin(i)
    lrootstot = lrootstot+lroots(i)*ispin(i)
  else
    lrootstot = lrootstot+lroots(i)
    nconftot = nconftot+nconf(i)
  end if
end do
! determine if number of roots equal to number of CSFs
if (lrootstot < nconftot) runmode = 4

! filling in lists of properties of states
call mma_allocate(list_sf_states,n_sf,label='list_sf_state')
call mma_allocate(list_sf_mult,n_sf,label='list_sf_mult')
call mma_allocate(list_sf_spin,n_sf,label='list_sf_spin')
call mma_allocate(list_so_mult,Nstate,label='list_so_mult')
call mma_allocate(list_so_spin,Nstate,label='list_so_spin')
call mma_allocate(list_so_proj,Nstate,label='list_so_proj')
call mma_allocate(list_so_sf,Nstate,label='list_so_sf')
ii = 1
jj = 1
do i=1,N ! spin manifolds
!sf values:
  do j=1,lroots(i)
    list_sf_states(ii) = j
    list_sf_mult(ii) = ispin(i)
    list_sf_spin(ii) = real((ispin(i)-One),kind=wp)/Two
    if (flag_so) then
      do m=1,ispin(i)
        list_so_mult(jj) = ispin(i)
        list_so_spin(jj) = list_sf_spin(ii)
        list_so_proj(jj) = real(m,kind=wp)-list_so_spin(jj)-1
        list_so_sf(jj) = ii
        jj = jj+1
      end do
    end if
    ii = ii+1
  end do
end do

if (ipglob > 2) then
  write(u6,*) 'sf_states: ',list_sf_states
  write(u6,*) 'sf_mult: ',list_sf_mult
  write(u6,*) 'sf_spin: ',list_sf_spin
  if (flag_so) then
    write(u6,*) 'so_mult: ',list_so_mult
    write(u6,*) 'so_proj: ',list_so_proj
    write(u6,*) 'so_sf: ',list_so_sf
    write(u6,*) 'so_spin: ',list_so_spin
  end if
end if

call mma_allocate(dipole,lrootstot,lrootstot,3,label='dipole')
call mma_allocate(dysamp,lrootstot,lrootstot,label='dysamp')
call mma_allocate(dysamp_bas,lrootstot,lrootstot,label='dysamp_bas')
call mma_allocate(tmp,Nstate,Nstate,label='tmp')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start from rassf/rassi output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if ((runmode /= 2) .and. (runmode /= 4)) then

  n_freq = Nstate*(Nstate-1)/2

  call mma_allocate(H_CSF,maxnconf,maxnconf,N)
  call mma_allocate(CI,maxnconf,maxlroots,N)
  call mma_allocate(E,maxlroots,N)
  call mma_allocate(U_CI,nconftot,lrootstot)
  call mma_allocate(DTOC,maxnum,maxnconf,N)
  call mma_allocate(V_CSF,nconftot,nconftot)
  call mma_allocate(HTOTRE_CSF,nconftot,nconftot)
  call mma_allocate(HTOT_CSF,nconftot,nconftot)
  call mma_allocate(a_einstein,lrootstot,lrootstot)
  call mma_allocate(emiss,n_freq)
  if (flag_so) then
    call mma_allocate(SO_CI,lrootstot,lrootstot)
    call mma_allocate(V_SO,lrootstot,lrootstot)
    call mma_allocate(CSF2SO,nconftot,lrootstot)
    call mma_allocate(E_SO,lrootstot)
  end if
  if (basis == 'SPH') call mma_allocate(E_SF,n_sf)

  ! Create PREP file for storage of intermediate data
  call cre_prep()

  ! reading wavefunction expansion from rasscf files
  H_CSF = Zero
  CI = Zero
  ! Determine file names
  ! Expected N rassd files and 1 rassisd file
  call mma_allocate(rassd_list,N)
  do i=1,N
    write(rassd_list(i),'(A5,I1)') 'RASSD',i
    if (ipglob > 2) write(u6,*) 'Reading ',rassd_list(i)
    call read_rassd(i)
  end do
  call mma_deallocate(rassd_list)

  if (i_rasscf == 1) then
    ! currently disabled
    ! Diagonalize H(RASSCF) and sort ascending the eigenvalue get the
    ! corresponding eigenvector (CI coefficients) to matrix CI(:,:,N)
    do i=1,N
      if (ipglob > 3) then
        write(u6,sint) 'H(RASSCF) in CSF basis of spin manifold:',i
        call dashes()
        do k=1,15
          write(u6,*) (H_CSF(k,j,i),j=1,15)
        end do
        call dashes()
      end if
      call mma_allocate(Ham,lroots(i),lroots(i))
      Ham(:,:) = H_CSF(1:nconf(i),1:nconf(i),i)
      call sortci(nconf(i),Ham,E(1:nconf(i),i),CI(1:nconf(i),1:lroots(i),i),ipglob)
      call mma_deallocate(Ham)
    end do
  else if ((i_rasscf == 2) .or. (i_rasscf == 3)) then
    ! Construct SF Hamiltonians from CIs and Es for each spin manifold
    do i=1,N
      do j=1,lroots(i)
        H_CSF(j,j,i) = E(j,i)
      end do
      ! H_CSF=CI*diag(E)*CI^T
      call mma_allocate(Ham,nconf(i),lroots(i))
      call mult(CI(1:nconf(i),1:lroots(i),i),H_CSF(1:lroots(i),1:lroots(i),i),Ham)
      call mult(Ham,CI(1:nconf(i),1:lroots(i),i),H_CSF(1:nconf(i),1:nconf(i),i),.false.,.true.)
      call mma_deallocate(Ham)
    end do
  end if

  call mh5_put_dset(prep_ci,CI)
  call mh5_put_dset(prep_hcsf,H_CSF)
  ! construct CI transformation matrix U_CI
  call uci()
  ! obtain pure SOC contribution to Hamiltonian
  call read_rassisd()
  if (flag_so) call get_vsoc()
  call get_hcsf()
  ! construct transformation matrices between bases
  if (flag_so) call soci()
  ! process dipole moment
  if (flag_so .and. basis /= 'SPH') call get_dipole()
  ! construct initial density matrix
  call get_dm0()

  ! close file PREP
  call mh5_close_file(prep_id)

  if (ipglob > 1) then
    call dashes()
    write(u6,*) 'Preparation finished successfully'
    call dashes()
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start from intermediate preparation file PREP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if (runmode == 2) then
  call mma_allocate(HTOT_CSF,nconftot,nconftot)
  call read_prep()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! case when number of states requested is less than number of CSFs
! Hamiltonian is read from RASSI and not constructed from RASSCF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if (runmode == 4) then
  call mma_allocate(HSOCX,Nstate,Nstate,label='HSOCX')
  call mma_allocate(CI,maxnconf,maxlroots,N,label='CI')
  call mma_allocate(U_CI,nconftot,lrootstot,label='U_CI')
  if (flag_so) then
    call mma_allocate(E_SO,Nstate,label='E_SO')
    call mma_allocate(SO_CI,lrootstot,lrootstot,label='SO_CI')
    call mma_allocate(CSF2SO,nconftot,lrootstot,label='CSF2SO')
  else
    call mma_allocate(E_SF,Nstate,label='E_SF')
  end if
  call get_dipole()
  if ((DM_basis == 'CSF_SO') .or. (DM_basis == 'SF') .or. (DM_basis == 'ALL') .or. (DM_basis == 'CSF_SF')) then
    CI = Zero
    ! Determine file names
    ! Expected N rassd files and 1 rassisd file
    call mma_allocate(rassd_list,N,label='rassd_list')
    do i=1,N
      write(rassd_list(i),'(A5,I1)') 'RASSD',i
      write(u6,*) rassd_list(i)
      ! c
      write(u6,*) 'entering read_rassd: ',i
      call read_rassd(i)
    end do
    call mma_deallocate(rassd_list)
    ! construct CI transformation matrix
    call uci()
    ! H_SOC hamiltonian from RASSI is read (in SF basis)
    call read_rassisd()
    ! construct the transformation matrix CSF2SO from SO to CSF
    if (flag_so) call mult(cmplx(U_CI,kind=wp),SO_CI,CSF2SO)
  else
    ! this is just caution condition to make sure that CM case
    ! was tested with given DM basis
    write(u6,*) 'WARNING!!! Take care of bases in CM case'
    call dashes()
    call abend()
  end if
  !call read_rassisd()
  call get_dm0()
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dynamics part starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (runmode /= 3) then
  ! determine preliminary dimension d of propagated matrices
  select case (basis)
    case ('CSF')
      d = nconftot
    case ('SF')
      d = lrootstot
    case ('SO')
      d = lrootstot
    case ('SPH')
      d = n_sf
  end select
  call mma_allocate(hamiltonian,d,d,label='hamiltonian')
  call mma_allocate(density0,d,d,label='density0')
  call mma_allocate(dipole_basis,d,d,3,label='dipole_basis')
  call mma_allocate(hamiltoniant,Nstate,Nstate,label='hamiltoniant')
  call mma_allocate(densityt,Nstate,Nstate,label='densityt')
  call mma_allocate(decay,Nstate,Nstate,label='decay')
  call mma_allocate(U_CI_compl,nconftot,lrootstot,label='U_CI_compl')

  if (runmode /= 4) then
    ! prepare density and hamiltonian in required basis to propagate with.
    ! here supposed that these matrices are prepared in CSF basis
    call hamdens()
  else
    !hamiltonian = HSOCX ! transform Hamiltonian to SO basis if requested
    if (flag_so .and. basis == 'SO') then
      call transform(HSOCX,SO_CI,hamiltonian)
    else
      hamiltonian(:,:) = HSOCX
    end if
    density0(:,:) = DM0
    dipole_basis(:,:,:) = dipole
    U_CI_compl(:,:) = cmplx(U_CI,kind=wp)
    if (flag_dyson) then
      do i=1,3
        dipole_basis(:,:,i) = dipole_basis(:,:,i)+alpha*dysamp_bas
      end do
    end if
  end if
  ! dynamics will be performed internally using matrices of dimension d
  if (basis == 'CSF') then
    d = nconftot
  else if ((basis == 'SF') .or. (basis == 'SO')) then
    if (Nstate == lrootstot) then
      d = lrootstot
    else if (Nstate < lrootstot) then
      ! not all states were requested for dynamics
      d = Nstate
      call cut_matrices()
    end if
  end if

  ! construct dissipation operator if requested
  if (flag_diss) then
    call mma_allocate(kab_basis,d,d)
    call mma_allocate(k_bar_basis,d,d)
    if (kext) then
      call k_external()
    else
      call kab()
    end if
    if (ipglob > 3) then
      lu = isFreeUnit(19)
      call molcas_open(lu,'Kab_basis_eV.dat')
      do i=1,d
        write(lu,*) (kab_basis(i,j)*autoev,j=1,d)
      end do
      close(lu)
      lu = isFreeUnit(20)
      call molcas_open(lu,'K_bar_basis_eV.dat')
      do i=1,d
        write(lu,*) (k_bar_basis(i,j)*autoev,j=1,d)
      end do
      close(lu)
    end if
  end if

  if (flag_decay .or. flag_dyson) call prepare_decay()

  if (basis /= 'SPH') then
    call propagate()
  else
    call propagate_sph()
  end if

end if

! put info for the test
if ((basis == 'SF')) then
  call mma_allocate(pop_sf,Nstate,label='pop_sf')
  pop_sf(:) = [(real(densityt(i,i)),i=1,Nstate)]
  call Add_Info('POP_SF',pop_sf,Nstate,3)
  call mma_deallocate(pop_sf)
end if

! closing and deallocation
call StatusLine('RhoDyn:','Close files and deallocate memory')

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

! allocated in prepare part
if (allocated(list_sf_states)) call mma_deallocate(list_sf_states)
if (allocated(list_sf_mult)) call mma_deallocate(list_sf_mult)
if (allocated(list_sf_spin)) call mma_deallocate(list_sf_spin)
if (allocated(list_so_mult)) call mma_deallocate(list_so_mult)
if (allocated(list_so_spin)) call mma_deallocate(list_so_spin)
if (allocated(list_so_proj)) call mma_deallocate(list_so_proj)
if (allocated(list_so_sf)) call mma_deallocate(list_so_sf)
if (allocated(V_CSF)) call mma_deallocate(V_CSF)
if (allocated(V_SO)) call mma_deallocate(V_SO)
if (allocated(H_CSF)) call mma_deallocate(H_CSF)
if (allocated(CI)) call mma_deallocate(CI)
if (allocated(DTOC)) call mma_deallocate(DTOC)
if (allocated(E)) call mma_deallocate(E)
if (allocated(U_CI)) call mma_deallocate(U_CI)
if (allocated(SO_CI)) call mma_deallocate(SO_CI)
if (allocated(E_SO)) call mma_deallocate(E_SO)
if (allocated(E_SF)) call mma_deallocate(E_SF)
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

call StatusLine('RhoDyn:','Finished')
ireturn = 0

return

end subroutine rhodyn
