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
!               2023, Thies Romig                                      *
!***********************************************************************

subroutine propagate_sph()
!***********************************************************************
! performs integration of LvN equation in the basis of irreducible
! tensors
!***********************************************************************

use integrators, only: rk4_sph
use rhodyn_data, only: d, densityt, dgl, dipole_basis, errorthreshold, finaltime, flag_fdm, flag_pulse, hamiltonian, hamiltoniant, &
                       hamiltoniant, initialtime, ipglob, ispin, k_max, k_ranks, len_sph, list_so_proj, list_so_sf, list_so_spin, &
                       list_sf_spin, lu_pls, lu_sf, method, midk1, midk2, midk3, midk4, mirr, n, nconftot, Nstate, Nstep, &
                       Ntime_tmp_dm, Npop, out3_fmt, out_fdmr, out_tfdm, q_max, q_proj, threshold, time_fdm, timestep, tout, v_so, &
                       v_so_red, Y1, Y2
use rhodyn_utils, only: check_hermicity, compare_matrices, dashes, print_c_matrix, werdm, werdm_back, werso, werso_back, W3J, W6J
use mh5, only: mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, cZero, auToFs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: a, b, c, i, ihh, ii, imm, iss, jj, k, K_prime, l, q, sa, sb, sc, m, Ntime, Noutstep
real(kind=wp) :: dum(3), time, timer(3), fact3j
character(len=64) :: sline
real(kind=wp), allocatable :: dgl_csf(:)
complex(kind=wp), allocatable :: density_csf(:,:), rho_init(:,:,:), rho_sph_t(:,:,:), tmp_back(:,:)
!procedure(pulse_func) :: pulse

call StatusLine('RhoDyn:','Propagation in Spherical Tensor basis starts')

if (ipglob > 2) call print_c_matrix(densityt,Nstate,'Initial density in SO basis')
! parameters of spherical decomposition
len_sph = (k_max+1)*(k_max+1)
if (q_max == -1) q_max = k_max
write(u6,'(a,i5,i5,i5,i5,a)') 'Dimension of the propagation: (',d,d,k_max+1,2*k_max+1,' )'
call mma_allocate(k_ranks,len_sph)
call mma_allocate(q_proj,len_sph)
call mma_allocate(rho_init,len_sph,d,d)
i = 1
do k=0,k_max
  do q=-k,k,1
    k_ranks(i) = k
    q_proj(i) = q
    ! initial density matrix decomposition
    call WERDM(densityt,Nstate,d,k,q,list_so_spin,list_so_proj,list_so_sf,rho_init(i,:,:))
    write(u6,'(a, i3, i3)') 'k, q: ',k,q
    if (ipglob > 2) call print_c_matrix(rho_init(i,:,:),d,'Density in ITO basis')
    i = i+1
  end do
end do
! check here back expansion of DM
call mma_allocate(tmp_back,Nstate,Nstate)
call WERDM_back(rho_init,Nstate,d,len_sph,k_ranks,q_proj,list_so_spin,list_so_proj,list_so_sf,tmp_back)
if (ipglob > 2) call print_c_matrix(tmp_back,Nstate,'Reexpanded initial density in SO basis')
call compare_matrices(densityt,tmp_back,Nstate,'Comparing DM decomposition',threshold)

! V_SO matrix decomposition
call mma_allocate(V_SO_red,d,d,3)
call WERSO(V_SO,Nstate,d,list_so_sf,list_so_spin,list_so_proj,V_SO_red)
if (ipglob > 2) then
  do m=1,3
    write(u6,*) 'm = ',m
    call print_c_matrix(V_SO_red(:,:,m),d,'Reduced V_SO')
  end do
end if
call WERSO_back(V_SO_red,Nstate,d,list_so_sf,list_so_spin,list_so_proj,tmp_back)
call compare_matrices(V_SO,tmp_back,Nstate,'Comparing V_SO decomposition',threshold)

if (ipglob > 2) then
  do m=1,3
    write(u6,*) 'm = ',m
    call print_c_matrix(dipole_basis(:,:,m),d,'Dipole')
  end do
end if

! prepare Y matrices
! structure should be the same as in equation_sph (last index is of <len_sph*3*n*2*(k_max+1) ??) (90 for ticl4)
call mma_allocate(Y1,d,d,1000)
call mma_allocate(Y2,d,d,1000)
Y1(:,:,:) = cZero
Y2(:,:,:) = cZero
i = 1
! l loop over indices k,q of ITOs basis
do l=1,len_sph
  k = k_ranks(l)
  q = q_proj(l)
  if ((q >= -q_max) .and. (q <= 0)) then
    do m=1,3
      do K_prime=k-1,k+1,1
        if ((K_prime < 0) .or. (K_prime > k_max)) cycle
        if ((q-m+2 > K_prime) .or. (q-m+2 < -K_prime)) cycle
        ! 3j symbol: (K_prime  1  k, q-m  m -q)
        fact3j = W3J(real(K_prime,kind=wp),One,real(k,kind=wp),real(q-m+2,kind=wp),real(m-2,kind=wp),real(-q,kind=wp))
        ! loop over spin manifolds
        do c=1,n
          sc = ispin(c)-1
          if (ipglob > 2) write(u6,*) 'i=',i
          ! a loop over rows
          do a=1,d
            sa = nint(2*list_sf_spin(a))
            ! b loop over columns
            do b=1,d
              sb = nint(2*list_sf_spin(b))
              !if (abs(V_SO_red(a,b,m)) < threshold) cycle
              Y1(a,b,i) = V_SO_red(a,b,m)*sqrt(real(3*(2*k+1)*(2*K_prime+1),kind=wp))*W6J(2*K_prime,2,2*k,sa,sc,sb)* &
                          fact3j*(-1)**((sa+sc)/2-q+m)
              Y2(a,b,i) = V_SO_red(a,b,m)*sqrt(real(3*(2*k+1)*(2*K_prime+1),kind=wp))*W6J(2,2*K_prime,2*k,sc,sb,sa)* &
                          fact3j*(-1)**((sb+sc)/2-q+m+K_prime+k)
            end do
          end do
          i = i+1
        end do
      end do
    end do
  end if
end do

! prepare mask matrix for mirroring
call mma_allocate(mirr,d,d)
do a=1,d
  sa = nint(2*list_sf_spin(a))
  do b=1,d
    sb = nint(2*list_sf_spin(b))
    mirr(a,b) = (-1)**((sa-sb)/2+1)
  end do
end do

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
  call mh5_put_dset(out_fdmr,abs(rho_init),[1,len_sph,d,d],[0,0,0,0])
end if

call mma_allocate(rho_sph_t,len_sph,d,d)
call mma_allocate(dgl,Nstate) ! accumulates diagonal elements of reexpanded matrix
call mma_allocate(midk1,len_sph,d,d)
call mma_allocate(midk2,len_sph,d,d)
call mma_allocate(midk3,len_sph,d,d)
call mma_allocate(midk4,len_sph,d,d)
if (allocated(hamiltoniant)) call mma_deallocate(hamiltoniant)
call mma_allocate(hamiltoniant,d,d,label='hamiltoniant')

rho_sph_t(:,:,:) = rho_init

! write initial values at 0th iteration
call mma_allocate(dgl_csf,nconftot,label='dgl_csf') ! dummy
call mma_allocate(density_csf,nconftot,nconftot,label='density_csf') ! dummy
call pop(time,ii,dgl_csf,density_csf)

do Ntime=1,(Nstep-1)
  write(sline,'(f10.3)') time*auToFs
  call StatusLine('RhoDyn: current time ',trim(sline))
  call Timing(dum(1),dum(2),timer(1),dum(3))
  if (flag_pulse) then
    ! update hamiltonian with dipole term
    call pulse(hamiltonian,hamiltoniant,time,Ntime)
  end if
  if (method == 'RK4_SPH') call rk4_sph(time,rho_sph_t)
  time = initialtime+timestep*Ntime
  if (mod(Ntime,Noutstep) == 0) then
    ii = ii+1
    ! transform density matrix back
    call WERDM_back(rho_sph_t,Nstate,d,len_sph,k_ranks,q_proj,list_so_spin,list_so_proj,list_so_sf,densityt)
    call pop(time,ii,dgl_csf,density_csf)
    if (flag_fdm .and. (time >= time_fdm*jj)) then
      ! write spherical density to file
      if (ipglob > 3) then
        i = 1
        do k=0,k_max
          do q=-k,k,1
            write(u6,*) 'time,k,q: ',time*auToFs,k,q
            call print_c_matrix(rho_sph_t(i,:,:),d,'Density in ITO basis')
            i = i+1
          end do
        end do
      end if
      call mh5_put_dset(out_tfdm,[time*auToFs],[1],[jj])
      call mh5_put_dset(out_fdmr,abs(rho_sph_t),[1,len_sph,d,d],[jj,0,0,0])
      jj = jj+1
    end if
  end if
  ! timing
  call Timing(dum(1),dum(2),timer(2),dum(3))
  timer(3) = timer(2)-timer(1)
  ihh = int(timer(3)/3600.0_wp)
  imm = int((timer(3)-ihh*3600.0_wp)/60.0_wp)
  iss = int(timer(3)-ihh*3600.0_wp-imm*60.0_wp)
  if (ipglob > 2) write(u6,out3_fmt) time*auToFs,time*auToFs,timestep*auToFs,ihh,':',imm,':',iss
  ! end timing
end do
call mma_deallocate(dgl_csf)
call mma_deallocate(density_csf)

call mma_deallocate(k_ranks)
call mma_deallocate(q_proj)
call mma_deallocate(rho_init)
call mma_deallocate(rho_sph_t)
call mma_deallocate(V_SO_red)
call mma_deallocate(tmp_back)
call mma_deallocate(midk1)
call mma_deallocate(midk2)
call mma_deallocate(midk3)
call mma_deallocate(midk4)
call mma_deallocate(dgl)
call mma_deallocate(Y1)
call mma_deallocate(Y2)
call mma_deallocate(mirr)

if (flag_pulse) close(lu_pls)
close(lu_sf)

call dashes()
write(u6,*) 'Propagation finished after ',Ntime,' steps'
call dashes()

call check_hermicity(densityt,Nstate,'Densityt',errorthreshold)

end subroutine propagate_sph
