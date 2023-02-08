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

subroutine propagate()
!***********************************************************************
! performs integration of LvN equation
!***********************************************************************

use integrators, only: classic_rk4, rk4, rk5, rk45, rkck
use rhodyn_data, only: ak1, ak2, ak3, ak4, ak5, ak6, d, decay, density0, densityt, dgl, DM_basis, dt, emiss, errorthreshold, &
                       finaltime, flag_decay, flag_dipole, flag_emiss, flag_fdm, flag_pulse, hamiltonian, hamiltoniant, &
                       initialtime, ipglob, lu_csf, lu_dip, lu_pls, lu_sf, lu_so, method, nconftot, Npop, Nstep, Ntime_tmp_dm, &
                       out2_fmt, out3_fmt, out_decay_i, out_decay_r, out_fdm, out_freq, out_ham_i, out_ham_r, out_tfdm, &
                       safety, time_fdm, timestep, tout
use rhodyn_utils, only: check_hermicity, dashes
use mh5, only: mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Five, Ten, Quart, auToFs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ihh, ii, imm, iss, jj, kk, Ntime, Noutstep
real(kind=wp) :: dum(3), error_rk, oldstep, t_temp, time, timer(3)
real(kind=wp), allocatable :: dgl_csf(:)
complex(kind=wp), allocatable :: density_csf(:,:)
!procedure(pulse_func) :: pulse

call dashes()
write(u6,*) 'Propagation starts'
write(u6,*) 'Dimension: ',d
call dashes()

! initialize parameters for solution of Liouville equation
ii = 1 ! counts output of populations
Ntime = 1 !counts steps
Nstep = int((finaltime-initialtime)/timestep)+1
Npop = int((finaltime-initialtime)/tout)+1
Noutstep = int(tout/timestep)
Ntime_tmp_dm = int(finaltime/time_fdm)+1 !fdm
time = initialtime
oldstep = timestep
densityt(:,:) = density0
! create and initialize h5 output file
call cre_out()
call mh5_put_dset(out_ham_r,real(hamiltonian))
call mh5_put_dset(out_ham_i,aimag(hamiltonian))
if (flag_decay) then
  call mh5_put_dset(out_decay_r,real(decay))
  call mh5_put_dset(out_decay_i,aimag(decay))
end if
if (flag_emiss) then
  call mh5_put_dset(out_freq,emiss)
  !write(u6,*) 'frequencies have been written to hdf5 file'
  emiss = Zero
end if
if (flag_fdm) then
  jj = 1 ! counts output of full density matrix
  ! store full density matrix
  call mh5_put_dset(out_tfdm,[time*auToFs],[1],[0])
  call mh5_put_dset(out_fdm,abs(density0),[1,d,d],[0,0,0])
end if

call mma_allocate(dgl,d,label='dgl')
call mma_allocate(ak1,d,d,label='ak1')
call mma_allocate(ak2,d,d,label='ak2')
call mma_allocate(ak3,d,d,label='ak3')
call mma_allocate(ak4,d,d,label='ak4')
call mma_allocate(ak5,d,d,label='ak5')
call mma_allocate(ak6,d,d,label='ak6')

call mma_allocate(dgl_csf,nconftot,label='dgl_csf')
call mma_allocate(density_csf,nconftot,nconftot,label='density_csf')
call pop(time,ii,dgl_csf,density_csf) ! write 0th iteration (initial values)

if ((method == 'RKCK') .or. (method == 'RK45')) then
  !*********************************************************************
  ! methods with adaptive step size
  !*********************************************************************
  kk = 1 ! counts initial steps by timestep
  ! initialize time step as 1/10 of interval:
  dt = (finaltime-initialtime)/Ten
  if (ipglob > 2) then
    call dashes()
    write(u6,*) 'Propagation with ',trim(method),' method'
    write(u6,*) 'Integration values at each step:'
    write(u6,out2_fmt) 'Time','Error','Timestep'
  end if
  ! main loop loopprop:
  loopprop: do while (time <= finaltime)
    call Timing(dum(1),dum(2),timer(1),dum(3))
    ! calculate hamiltonian with pulse at the current time:
    if (flag_pulse) then
      if (time >= (initialtime+timestep*(kk-1))) then
        call pulse(hamiltonian,hamiltoniant,time,kk)
        kk = kk+1
      else
        call pulse(hamiltonian,hamiltoniant,time,-1)
      end if
    else
      hamiltoniant(:,:) = hamiltonian
    end if
    ! run loop to find the step with acceptable accuracy
    ! use density0 as storage for initial value
    loopstep: do
      density0(:,:) = densityt
      if (method == 'RKCK') then
        call rkck(time,densityt,error_rk)
      else if (method == 'RK45') then
        call rk45(time,densityt,error_rk)
      end if
      if (error_rk <= errorthreshold) exit loopstep
      ! then step rejected, try new smaller step
      densityt(:,:) = density0
      t_temp = safety*dt*(errorthreshold/error_rk)**Quart
      dt = max(t_temp,0.2_wp*dt)
      if (ipglob > 2) write(u6,*) ' ------',error_rk,dt*auToFs
    end do loopstep
    ! step succeeded
    time = time+dt
    oldstep = dt
    if (time >= finaltime) exit loopprop
    ! write info and elapsed time
    if (time >= (initialtime+tout*ii)) then
      ii = ii+1
      call pop(time,ii,dgl_csf,density_csf)
    end if
    if (flag_fdm .and. (time >= time_fdm*jj)) then
      ! should be moved to procedure pop
      call mh5_put_dset(out_tfdm,[time*auToFs],[1],[jj])
      ! density0 is stored as temporary storage for dm in required basis in pop
      call mh5_put_dset(out_fdm,abs(density0),[1,d,d],[jj,0,0])
      jj = jj+1
    end if
    call Timing(dum(1),dum(2),timer(2),dum(3))
    timer(3) = timer(2)-timer(1)
    ihh = int(timer(3)/3600.0_wp)
    imm = int((timer(3)-ihh*3600.0_wp)/60.0_wp)
    iss = int(timer(3)-ihh*3600.0_wp-imm*60.0_wp)
    if (ipglob > 2) write(u6,out3_fmt) time*auToFs,error_rk,dt*auToFs,ihh,':',imm,':',iss
    ! next step size
    if (method == 'RKCK') then
      dt = safety*dt*(errorthreshold/error_rk)**0.2_wp
    else if (method == 'RK45') then
      dt = safety*dt*(errorthreshold/error_rk)**Quart
    end if
    ! stepsize cannot decrease faster than by 5x
    dt = min(dt,Five*oldstep)
    Ntime = Ntime+1
    if (dt < 1.0e-8_wp) then
      write(u6,*) 'The step of integration got too small, try fixed-step methods and/or check the hamiltonian'
      call abend()
    end if
  end do loopprop
else
  !*********************************************************************
  ! methods with fixed step
  !*********************************************************************
  do Ntime=1,(Nstep-1)
    if (flag_pulse) then
      ! update hamiltonian with dipole term
      call pulse(hamiltonian,hamiltoniant,time,Ntime)
    else
      hamiltoniant(:,:) = hamiltonian
    end if
    select case (method)
      case ('CLASSIC_RK4')
        call classic_rk4(time,densityt)
      case ('RK4')
        call rk4(time,densityt)
      case ('RK5')
        call rk5(time,densityt)
    end select
    time = initialtime+timestep*Ntime
    if (mod(Ntime,Noutstep) == 0) then
      ii = ii+1
      call pop(time,ii,dgl_csf,density_csf)
    end if
  end do
end if
call mma_deallocate(dgl_csf)
call mma_deallocate(density_csf)

! deallocation of matrices needed for propagation
call mma_deallocate(ak1)
call mma_deallocate(ak2)
call mma_deallocate(ak3)
call mma_deallocate(ak4)
call mma_deallocate(ak5)
call mma_deallocate(ak6)
call mma_deallocate(dgl)

if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) close(lu_so)
if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) close(lu_sf)
if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) close(lu_csf)
if (flag_pulse) close(lu_pls)
if (flag_dipole) close(lu_dip)

call dashes()
write(u6,*) 'Propagation finished after ',Ntime,' steps'
call dashes()

call check_hermicity(densityt,d,'Densityt',errorthreshold)

end subroutine propagate
