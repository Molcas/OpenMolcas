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

use Constants, only: Zero, One, cZero, cOne, auToFs
use Definitions, only: wp, iwp, u6
use integrators, only: rk4_sph
use rhodyn_data, only: d, dgl, dipole_basis, E_SF, finaltime, flag_pulse, hamiltonian, initialtime, ipglob, len_sph, list_so_proj, &
                       list_so_sf, list_so_spin, method, midk1, midk2, midk3, midk4, n_sf, N_Populated, Nstate, Nstep, k_max, &
                       k_ranks, pulse_func, q_proj, rho_sph_t, threshold, timestep, tout, V_SO, V_SO_red
use rhodyn_utils, only: dashes, werdm, WERDM_back, WERSO, WERSO_back, print_c_matrix, compare_matrices, check_hermicity
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer(kind=iwp) :: i, k, q, m, Ntime, Noutstep, lu
real(kind=wp) :: time
complex(kind=wp), allocatable :: dm0_so(:,:), dm0_so_back(:,:), V_SO_back(:,:)
complex(kind=wp), allocatable :: rho_init(:,:,:), dum_zero(:,:)
integer(kind=iwp), external :: isFreeUnit
procedure(pulse_func) :: pulse
character(len=64) :: sline, out_fmt_sph

! so density matrix preparation
call mma_allocate(dm0_so,Nstate,Nstate)
call mma_allocate(dm0_so_back,Nstate,Nstate)
dm0_so = cZero
dm0_so(N_Populated,N_Populated) = cOne
if (ipglob>2) call print_c_matrix(dm0_so,Nstate,'Initial density in SO basis')
! parameters of spherical decomposition
! k_max should be defined
!k_max = 2
len_sph = (k_max+1)*(k_max+1)
!call mma_allocate(k_ranks,k_max+1)
call mma_allocate(k_ranks,len_sph)
call mma_allocate(q_proj,len_sph)
call mma_allocate(rho_init,len_sph,n_sf,n_sf)
i = 1
do k=0,k_max
  !k_ranks(k+1) = k
  do q=-k,k,1
    k_ranks(i)= k
    q_proj(i) = q
    ! initial density matrix decomposition
    call WERDM(dm0_so,Nstate,n_sf,k,q,list_so_spin,list_so_proj,list_so_sf,rho_init(i,:,:))
    write(u6,*) 'time,k,q: ', Zero, k, q
    if (ipglob>2) call print_c_matrix(rho_init(i,:,:), n_sf, 'Density in ITO basis')
    i=i+1
  end do
end do
! check here back expansion of DM
call WERDM_back(rho_init,Nstate,n_sf,len_sph,k_ranks,q_proj,list_so_spin,list_so_proj,list_so_sf,dm0_so_back)
if (ipglob > 2) call print_c_matrix(dm0_so_back,Nstate,'Reexpanded initial density in SO basis')
call compare_matrices(dm0_so, dm0_so_back, Nstate, 'Comparing DM decomposition', threshold)

write(u6,*) 'k_ranks: ', k_ranks
write(u6,*) 'q_proj: ', q_proj

! working with V_SO
call mma_allocate(V_SO_red,n_sf,n_sf,3)
call mma_allocate(V_SO_back,Nstate,Nstate)
call WERSO(V_SO,Nstate,n_sf,list_so_sf,list_so_spin,list_so_proj,V_SO_red)
if (ipglob > 2) then
  do m=1,3
    write(u6,*) 'm = ', m
    call print_c_matrix(V_SO_red(:,:,m), n_sf, 'Reduced V_SO')
  end do
end if
call WERSO_back(V_SO_red,Nstate,n_sf,list_so_sf,list_so_spin,list_so_proj,V_SO_back)
call compare_matrices(V_SO, V_SO_back, Nstate, 'Comparing V_SO decomposition', threshold)

if (ipglob > 2) then
  write(u6,*) 'Energies: ', E_SF
  do m=1,3
    write(u6,*) 'm = ', m
    call print_c_matrix(dipole_basis(:,:,m), n_sf, 'Dipole')
  end do
end if

call dashes()
write(u6,*) 'Propagation starts'
write(u6,*) 'Dimension: (', d, d, k_max, 2*k_max+1, ' )'
call dashes()

! initialize parameters for solution of Liouville equation
Nstep = int((finaltime-initialtime)/timestep)+1
Noutstep = int(tout/timestep)
time = initialtime

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

lu = isFreeUnit(11)
call molcas_open(lu,'density_so.out')
dgl(:) = [(real(dm0_so_back(i,i)),i=1,Nstate)]
write(out_fmt_sph,'(a,i5,a)') '(1x,',Nstate+1,'(f22.16))'
write(lu,out_fmt_sph) time*auToFs,(dgl(i),i=1,Nstate)

do Ntime=1,(Nstep-1)
  write(sline,'(f10.3)') time*auToFs
  call StatusLine('RhoDyn: current time ',trim(sline))
  if (flag_pulse) then
    ! update hamiltonian with dipole term
    call pulse(dum_zero,hamiltonian,time,-1)
  end if
  if (method == 'RK4_SPH') call rk4_sph(time,rho_sph_t)
!   !!vk!! call test_rho(densityt,time)
  time = initialtime+timestep*Ntime
  if (mod(Ntime,Noutstep) == 0) then
    !transform matrix back, using dm0_so_back as temp storage
    call WERDM_back(rho_sph_t,Nstate,n_sf,len_sph,k_ranks,q_proj,list_so_spin,list_so_proj,list_so_sf,dm0_so_back)
    dgl(:) = [(real(dm0_so_back(i,i)),i=1,Nstate)]
    write(lu,out_fmt_sph) time*auToFs,(dgl(i),i=1,Nstate)
    if (ipglob > 3) then
      i = 1
      do k=0,k_max
        do q=-k,k,1
          write(u6,*) 'time,k,q: ', time*auToFs, k, q
          call print_c_matrix(rho_sph_t(i,:,:), n_sf, 'Density in ITO basis')
          i=i+1
        end do
      end do
    end if
  end if
end do

! closing files and deallocation
close(lu)

if (allocated(k_ranks)) call mma_deallocate(k_ranks)
if (allocated(q_proj)) call mma_deallocate(q_proj)
if (allocated(dm0_so)) call mma_deallocate(dm0_so)
if (allocated(dm0_so_back)) call mma_deallocate(dm0_so_back)
if (allocated(rho_init)) call mma_deallocate(rho_init)
if (allocated(V_SO_red)) call mma_deallocate(V_SO_red)
if (allocated(V_SO_back)) call mma_deallocate(V_SO_back)
if (allocated(dum_zero)) call mma_deallocate(dum_zero)
call mma_deallocate(midk1)
call mma_deallocate(midk2)
call mma_deallocate(midk3)
call mma_deallocate(midk4)
call mma_deallocate(dgl)

end subroutine propagate_sph
