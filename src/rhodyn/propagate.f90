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
subroutine propagate()
  use rhodyn_data
  use rhodyn_utils, only: dashes
  use definitions, only: wp, iwp, u6
  use constants, only: auToFs
  use stdalloc, only: mma_allocate, mma_deallocate
  use mh5, only: mh5_put_dset
  implicit none
!***********************************************************************
!
!
!
!***********************************************************************
  integer(kind=iwp) :: Ntime,Noutstep
! timers are commented out now, should be again switched on, once it is
! transferred to f90 standard
!  integer :: ihh, imm, iss
  real(kind=wp) :: time, error_rk, oldstep, t_temp
  integer(kind=iwp), external :: isfreeunit
  procedure(rk_fixed_step) :: classic_rk4,rk4,rk5
  procedure(rk_adapt_step) :: rk45,rkck
  procedure(pulse_func)    :: pulse
  character(len=64) :: out2_fmt, out3_fmt

! timers are commented out now, should be again switched on, once it is
! transferred to f90 standard
!#include "timers.fh"

  call dashes()
  write(u6,*) 'Propagation starts'
  call dashes()

! write formats for output files SFDENS, SODENS, CSFDENS
  ! header format
  write(out1_fmt,"(a,i5,a)") "(a22,",Nstate,"(i8,14x),a)"
  write(out1_fmt_csf,"(a,i5,a)") "(a22,",nconftot,"(i8,14x),a)"
  ! line format
  write(out_fmt, "(a,i5,a)") "(x,",Nstate+2,"(f22.16))"
  write(out_fmt_csf, "(a,i5,a)") "(x,",nconftot+2,"(f22.16))"
  out2_fmt='(2x,a,28x,a,28x,a,28x,a)'
  out3_fmt='(2x,f22.16,1x,f22.16,1x,f22.16,1x,i1,a1,i2.2,a1,i2.2)'

! this file is for the diagonal density matrix in SO basis
  if ((DM_basis=='SO').or.(DM_basis=='CSF_SO').or. &
      (DM_basis=='SF_SO').or.(DM_basis=='ALL')) then
    lu_so=isfreeunit(lu_so)
    call molcas_open (lu_so,'SODENS')
    write(lu_so,out1_fmt) '#time(fs)',(j,j=1,Nstate),'Norm'
  endif
! this file is for the diagonal density matrix in SF basis
  if ((DM_basis=='SF').or.(DM_basis=='CSF_SF').or. &
      (DM_basis=='SF_SO').or.(DM_basis=='ALL')) then
    lu_sf=isfreeunit(lu_sf)
    call molcas_open (lu_sf,'SFDENS')
    write(lu_sf,out1_fmt) '#time(fs)',(j,j=1,Nstate),'Norm'
  endif
! this file is for the diagonal density matrix in CSF basis
  if ((DM_basis=='CSF').or.(DM_basis=='CSF_SF').or. &
      (DM_basis=='CSF_SO').or.(DM_basis=='ALL')) then
    lu_csf=isfreeunit(lu_csf)
    call molcas_open (lu_csf,'CSFDEN')
    write(lu_csf,out1_fmt_csf) '#time(fs)',(j,j=1,nconftot),'Norm'
  endif
! this file is for the pulse data
  if (flag_pulse) then
    lu_pls=isfreeunit(lu_pls)
    call molcas_open(lu_pls,'PULSE')
  endif
! this file is for the TD-dipole moment data
  if (flag_dipole) then
    lu_dip=isfreeunit(lu_dip)
    call molcas_open(lu_dip,'DIPOLE.dat')
  endif

! initialize parameters for solution of Liouville equation
  Nstep   = int((finaltime-initialtime)/timestep)+1
  Npop    = int((finaltime-initialtime)/tout)+1
  Noutstep= int(tout/timestep)
  if (flag_fdm) Ntime_tmp_dm=int(finaltime/time_fdm)+1
  time    = initialtime
  oldstep = timestep
  densityt= density0
  call cre_out()
  call mh5_put_dset(out_ham_r,dble(hamiltonian))
  call mh5_put_dset(out_ham_i,aimag(hamiltonian))
  call mh5_put_dset(out_decay_r,dble(decay))
  call mh5_put_dset(out_decay_i,aimag(decay))
  if (flag_emiss) then
    call mh5_put_dset(out_freq, emiss)
!   write(*,*) 'frequencies have been written to hdf5 file'
    emiss = 0d0
  endif
  ii = 1 ! counts output of populations
  jj = 1 ! counts output of full density matrix (flag_fdm)
  Ntime=1 !counts steps
  call mma_allocate(dgl,d)
  call pop(time,ii)
  if (flag_fdm) then
    ! store full density matrix
    call mh5_put_dset(out_tfdm, [time*auToFs], [1], [0])
    call mh5_put_dset(out_fdm,abs(density0),[1,d,d],[0,0,0])
  endif
  call mma_allocate(ak1,d,d)
  call mma_allocate(ak2,d,d)
  call mma_allocate(ak3,d,d)
  call mma_allocate(ak4,d,d)
  call mma_allocate(ak5,d,d)
  call mma_allocate(ak6,d,d)

!***********************************************************************
!     methods with adaptive step size
!***********************************************************************
  if (method=='RKCK'.or.method=='RK45') then
    kk=1 ! counts initial steps by timestep
    ! initialize time step as 1/10 of interval:
    dt=(finaltime-initialtime)/10.
    if (ipglob>2) then
      call dashes()
      write(u6,*) 'Propagation with ',trim(method),' method'
      write(u6,*) 'Integration values at each step:'
      write(u6,out2_fmt) 'Time','Error','Timestep'
    endif
    ! main loop loopprop:
    loopprop: do while (time<=finaltime)
!      call Timing(Swatch,Swatch,Rolex_1,Swatch)
!      calculate hamiltonian with pulse at the current time:
      if (flag_pulse) then
        if (time>=(initialtime+timestep*(kk-1))) then
          call pulse(hamiltonian,hamiltoniant,time,kk)
          kk=kk+1
        else
          call pulse(hamiltonian,hamiltoniant,time)
        endif
      else
        hamiltoniant=hamiltonian
      endif
! run loop to find the step with acceptable accuracy
! use density0 as storage for initial value
      loopstep: do
!       call dcopy_(d**2,densityt,1,density0,1)
        density0=densityt
        if (method=='RKCK') then
          call rkck(time,densityt,error_rk)
        elseif (method=='RK45') then
          call rk45(time,densityt,error_rk)
        endif
        if (error_rk<=errorthreshold) exit loopstep
! then step rejected, try new smaller step
        densityt=density0
!        call dcopy_(d**2,density0,1,densityt,1)
        t_temp=safety*dt*(errorthreshold/error_rk)**0.25
        dt=max(t_temp,0.2*dt)
        if (ipglob>2) write(u6,*) ' ------',error_rk,dt*auToFs
      enddo loopstep
! step succeeded
      time=time+dt
      oldstep=dt
      if (time>=finaltime) exit loopprop
! write info and elapsed time
      if (time>=(initialtime+tout*ii)) then
        ii=ii+1
        call pop(time,ii)
      endif
      if (flag_fdm.and.time>=time_fdm*jj) then
        ! should be moved to procedure pop
        call mh5_put_dset(out_tfdm,[time*auToFs],[1],[jj])
        call mh5_put_dset(out_fdm,abs(densityt), &
                                  [1,d,d],[jj,0,0])
            jj = jj + 1
      endif
!      call Timing(Swatch,Swatch,Rolex_2,Swatch)
!      Rolex_3=Rolex_2-Rolex_1
!      ihh=int(Rolex_3/3600)
!      imm=int(Rolex_3-ihh*3600)/60
!      iss=int(Rolex_3-ihh*3600-imm*60)
      if (ipglob>2) write(u6,out3_fmt) time*auToFs,error_rk,dt*auToFs
!                            , ihh,':',imm,':',iss
!      next step size
      if (method=='RKCK') then
        dt=safety*dt*(errorthreshold/error_rk)**0.2
      elseif (method=='RK45') then
        dt=safety*dt*(errorthreshold/error_rk)**0.25
      endif
      dt=min(dt,5*oldstep)
      Ntime=Ntime+1
      if (dt<1d-8) then
        write(u6,*) "The step of integration got too small," &
               //"try fixed-step methods and/or check the hamiltonian"
        call abend()
      endif
    enddo loopprop
  else
!***********************************************************************
!       methods with fixed step
!***********************************************************************
    do Ntime=1,(Nstep-1)
      if (flag_pulse) then
        call pulse(hamiltonian,hamiltoniant,time,Ntime)
      else
        hamiltoniant=hamiltonian
      endif
      select case (method)
      case ('CLASSIC_RK4')
        call classic_rk4(time,densityt)
      case ('RK4')
        call rk4(time,densityt)
      case ('RK5')
        call rk5(time,densityt)
      case default
        ! check has already been done in read_input.f90
        write(6,*)'Integration method ',method, ' is not included'
        call abend()
      end select
!!vk!!    call test_rho(densityt,time)
      time=initialtime+timestep*Ntime
      if (mod(Ntime,Noutstep)==0) then
        ii = ii + 1
        call pop(time,ii)
      endif
    enddo
  endif

  call mma_deallocate(ak1)
  call mma_deallocate(ak2)
  call mma_deallocate(ak3)
  call mma_deallocate(ak4)
  call mma_deallocate(ak5)
  call mma_deallocate(ak6)
  call mma_deallocate(dgl)

  if ((DM_basis=='SO').or.(DM_basis=='CSF_SO').or. &
         (DM_basis=='SF_SO').or.(DM_basis=='ALL')) close(lu_so)
  if ((DM_basis=='SF').or.(DM_basis=='CSF_SF').or. &
         (DM_basis=='SF_SO').or.(DM_basis=='ALL')) close(lu_sf)
  if ((DM_basis=='CSF').or.(DM_basis=='CSF_SF').or. &
         (DM_basis=='CSF_SO').or.(DM_basis=='ALL')) close(lu_csf)
  if (flag_pulse) close(lu_pls)
  if (flag_dipole) close(lu_dip)


  call dashes()
  write(6,*) 'Propagation finished after ', Ntime, ' steps'
  call dashes()

end
