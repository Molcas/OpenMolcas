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
! Copyright (C) 2022-2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine propagate_sph()
!***********************************************************************
! performs integration of LvN equation in the basis of irreducible
! tensors
!***********************************************************************

use integrators, only: rk4_sph
use rhodyn_data, only: d, densityt, dgl, dipole_basis, E_SF, errorthreshold, finaltime, flag_fdm, flag_pulse, &
                       hamiltonian, initialtime, ipglob, len_sph, list_so_proj, list_so_sf, list_so_spin, method, midk1, &
                       midk2, midk3, midk4, n_sf, N_Populated, nconftot, Nstate, Nstep, Ntime_tmp_dm, Npop, k_max, k_ranks, &
                       out_fdm, out_tfdm, q_proj, threshold, time_fdm, timestep, tout, V_SO, V_SO_red

use rhodyn_utils, only: dashes, werdm, WERDM_back, WERSO, WERSO_back, print_c_matrix, check_hermicity, compare_matrices
use mh5, only: mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne, auToFs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, ii, jj, k, q, m, Ntime, Noutstep
real(kind=wp) :: time
character(len=64) :: sline
real(kind=wp), allocatable :: dgl_csf(:)
complex(kind=wp), allocatable :: density_csf(:,:), dum_zero(:,:), rho_init(:,:,:), rho_sph_t(:,:,:), tmp_back(:,:)
!procedure(pulse_func) :: pulse

call StatusLine('RhoDyn:','Propagation in Spherical Tensor basis starts')

! so density matrix preparation
densityt = cZero
densityt(N_Populated,N_Populated) = cOne
if (ipglob > 2) call print_c_matrix(densityt,Nstate,'Initial density in SO basis')
! parameters of spherical decomposition
len_sph = (k_max+1)*(k_max+1)
write(u6,'(a,i5,i5,i5,i5,a)') 'Dimension of the propagation: (',d,d,k_max+1,2*k_max+1,' )'
call mma_allocate(k_ranks,len_sph)
call mma_allocate(q_proj,len_sph)
call mma_allocate(rho_init,len_sph,n_sf,n_sf)
i = 1
do k=0,k_max
  do q=-k,k,1
    k_ranks(i) = k
    q_proj(i) = q
    ! initial density matrix decomposition
    call WERDM(densityt,Nstate,n_sf,k,q,list_so_spin,list_so_proj,list_so_sf,rho_init(i,:,:))
    write(u6,'(a, i3, i3)') 'k, q: ',k,q
    if (ipglob > 2) call print_c_matrix(rho_init(i,:,:),n_sf,'Density in ITO basis')
    i = i+1
  end do
end do
! check here back expansion of DM
call mma_allocate(tmp_back,Nstate,Nstate)
call WERDM_back(rho_init,Nstate,n_sf,len_sph,k_ranks,q_proj,list_so_spin,list_so_proj,list_so_sf,tmp_back)
if (ipglob > 2) call print_c_matrix(tmp_back,Nstate,'Reexpanded initial density in SO basis')
call compare_matrices(densityt,tmp_back,Nstate,'Comparing DM decomposition',threshold)

! V_SO matrix decomposition
call mma_allocate(V_SO_red,n_sf,n_sf,3)
call WERSO(V_SO,Nstate,n_sf,list_so_sf,list_so_spin,list_so_proj,V_SO_red)
if (ipglob > 2) then
  do m=1,3
    write(u6,*) 'm = ',m
    call print_c_matrix(V_SO_red(:,:,m),n_sf,'Reduced V_SO')
  end do
end if
call WERSO_back(V_SO_red,Nstate,n_sf,list_so_sf,list_so_spin,list_so_proj,tmp_back)
call compare_matrices(V_SO,tmp_back,Nstate,'Comparing V_SO decomposition',threshold)

if (ipglob > 2) then
  write(u6,*) 'Energies: ',E_SF
  do m=1,3
    write(u6,*) 'm = ',m
    call print_c_matrix(dipole_basis(:,:,m),n_sf,'Dipole')
  end do
end if

! initialize parameters for solution of Liouville equation
ii = 1 ! counts output of populations
Nstep = int((finaltime-initialtime)/timestep)+1
Npop = int((finaltime-initialtime)/tout)+1
Noutstep = int(tout/timestep)
Ntime_tmp_dm = int(finaltime/time_fdm)+1 !fdm
time = initialtime
! create and initialize h5 output file
call cre_out()
if (flag_fdm) then
  jj = 1 ! counts output of full density matrix
  ! store full density matrix
  call mh5_put_dset(out_tfdm,[time*auToFs],[1],[0])
  call mh5_put_dset(out_fdm,abs(rho_init),[1,len_sph,d,d],[0,0,0,0])
end if

call mma_allocate(rho_sph_t,len_sph,n_sf,n_sf)
call mma_allocate(dgl,Nstate) ! accumulates diagonal elements of reexpanded matrix
call mma_allocate(midk1,len_sph,d,d)
call mma_allocate(midk2,len_sph,d,d)
call mma_allocate(midk3,len_sph,d,d)
call mma_allocate(midk4,len_sph,d,d)
if (flag_pulse) then
  call mma_allocate(dum_zero,d,d,label='dum_zero')
  dum_zero(:,:) = cZero
end if

rho_sph_t(:,:,:) = rho_init

! write initial values at 0th iteration
call mma_allocate(dgl_csf,nconftot,label='dgl_csf') ! dummy
call mma_allocate(density_csf,nconftot,nconftot,label='density_csf') ! dummy
call pop(time,ii,dgl_csf,density_csf)

do Ntime=1,(Nstep-1)
  write(sline,'(f10.3)') time*auToFs
  call StatusLine('RhoDyn: current time ',trim(sline))
  if (flag_pulse) then
    ! update hamiltonian with dipole term
    call pulse(dum_zero,hamiltonian,time,Ntime)
  end if
  if (method == 'RK4_SPH') call rk4_sph(time,rho_sph_t)
  time = initialtime+timestep*Ntime
  if (mod(Ntime,Noutstep) == 0) then
    ii = ii+1
    ! transform density matrix back
    call WERDM_back(rho_sph_t,Nstate,n_sf,len_sph,k_ranks,q_proj,list_so_spin,list_so_proj,list_so_sf,densityt)
    call pop(time,ii,dgl_csf,density_csf)
    if (flag_fdm .and. (time >= time_fdm*jj)) then
      ! write spherical density to file
      if (ipglob > 3) then
        i = 1
        do k=0,k_max
          do q=-k,k,1
            write(u6,*) 'time,k,q: ',time*auToFs,k,q
            call print_c_matrix(rho_sph_t(i,:,:),n_sf,'Density in ITO basis')
            i = i+1
          end do
        end do
      end if
      call mh5_put_dset(out_tfdm,[time*auToFs],[1],[jj])
      call mh5_put_dset(out_fdm,abs(rho_sph_t),[1,len_sph,d,d],[jj,0,0,0])
      jj = jj+1
    end if
  end if
end do
call mma_deallocate(dgl_csf)
call mma_deallocate(density_csf)

call mma_deallocate(k_ranks)
call mma_deallocate(q_proj)
call mma_deallocate(rho_init)
call mma_deallocate(rho_sph_t)
call mma_deallocate(V_SO_red)
call mma_deallocate(tmp_back)
if (allocated(dum_zero)) call mma_deallocate(dum_zero)
call mma_deallocate(midk1)
call mma_deallocate(midk2)
call mma_deallocate(midk3)
call mma_deallocate(midk4)
call mma_deallocate(dgl)

call dashes()
write(u6,*) 'Propagation finished after ',Ntime,' steps'
call dashes()

call check_hermicity(densityt,Nstate,'Densityt',errorthreshold)

end subroutine propagate_sph
