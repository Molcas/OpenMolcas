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
subroutine read_input()
  use rhodyn_data
  use rhodyn_utils, only: dashes
  use definitions, only: iwp, u6
  use constants, only: auToFs, auToCm, auToeV, pi
  use stdalloc, only: mma_allocate, mma_deallocate
  implicit none
!
!
!
  integer(kind=iwp)    :: luin
  character(len=8)     :: input_id = '&RHODYN'
  character(len=32)    :: tryname
  character(len=256)   :: line

!  call qEnter('read_input') ! deprecated

  call SpoolInp(luin)
! Find beginning of input:
50 read(luin,'(A72)') line
  call normal(line)
  if(line(1:8)/=input_id) goto 50

  do
    read(luin,'(A72)',end=300) line
    call normal(line)
    if(line(1:1)=='*') cycle
    if(line==' ') cycle
    select case (line(1:4))
    case('NRSM')
      read(luin,*)N
      allocate(ndet(N))
      allocate(nconf(N))
      allocate(lroots(N))
      allocate(ispin(N))
      case('NRDE')
        do i=1,N
          read(luin,*) ndet(i),nconf(i),lroots(i),ispin(i)
        enddo
      case('POPU')
        read(luin,'(A)') p_style
        call upCase(p_style)
        if (p_style/='DET'.and.p_style/='CSF'.and.p_style/='SF'.and.&
            p_style/='SO'.and.p_style/='SO_THERMAL'.and.&
            p_style/='SF_THERMAL'.and.p_style/='FROMFILE') then
          call WarningMessage(2,'Unknown option for POPUlation style')
          call abend()
        endif
      case('NRPO')
        read(luin,*) N_Populated
      case('TEMP')
        read(luin,*) T
      case('IFSO')
        flag_so=.True.
      case('NMOD')
        read(luin,'(I8)') Nmode
      case('PROP')
        read(luin,'(A)') basis
        call upCase(basis)
        if (basis/='DET'.and.basis/='CSF'.and.basis/='SF'.and.&
            basis/='SO') then
          call WarningMessage(2,'Unknown option for PROPagation basis')
          call abend()
        endif
      case('NSTA')
        read(luin,*) Nstate,tryname
        call UpCase(tryname)
        if (tryname=='ALL') then
          istates=(/(i,i=1,Nstate)/)
        else
          backspace(luin)
          allocate(istates(Nstate))
          read(luin,*) Nstate, (istates(i),i=1,Nstate)
        endif
      case('PREP')
        read(luin,'(I8)') preparation
      case('TOUT')
        read(luin,*) tout
        tout=tout/auToFs
      case('INIT')
        read(luin,*) initialtime
        initialtime=initialtime/auToFs
      case('FINA')
        read (luin,*) finaltime
        finaltime=finaltime/auToFs
      case('TSTE')
        read(luin,*) timestep
        timestep=timestep/auToFs
      case('METH')
        read(luin,*) method
        call upCase(method)
        if (method/='RKCK'.and.method/='RK45'.and.&
            method/='RK4'.and.method/='RK5'.and.&
            method/='CLASSIC_RK4') then
          call WarningMessage(2,'Unknown option for METHod')
          call abend()
        endif
      case('RK45')
        read(luin,*) errorthreshold
      case('RKSA')
        read(luin,*) safety
      case('DELT')
        read(luin,*) deltaE
        deltaE=deltaE/auToCm
      case ('VCOU')
        read (luin,*) V
        V=V/auToCm
!       case ('Out_fmt')
!         read(11,'(A)') Out_fmt
      case ('AUGE')
        flag_decay=.True.
      case ('NVAL')
        read(luin,'(I8)') Nval
      case ('DECA')
        read(luin,*) N_L3, tau_L3
        read(luin,*) N_L2, tau_L2
        tau_L3=tau_L3/autoev
        tau_L2=tau_L2/autoev
!     tau_L3=tau_L3/auToFs
!     tau_L2=tau_L2/auToFs
      case ('DYSO')
        flag_dyson=.True.
      case ('ALPH')
        read(luin,*) alpha
      case ('IFDI')
        flag_diss=.True.
      case ('IOND')
        read(luin,*) ion_diss
      case ('GAMM')
        read(luin,*) gamma
        gamma=gamma/auToCm
      case ('HRSO')
        HRSO=.True.
      case ('KEXT')
        kext=.True.
      case ('TFDM')
        read(luin,*) time_fdm
        time_fdm=time_fdm/auToFs
        flag_fdm = .True.
      case ('DMBA')
        read(luin,'(A)') dm_basis
        call upCase(dm_basis)
        if (dm_basis/='CSF'.and.dm_basis/='SF'.and.dm_basis/='SO'.and.&
            dm_basis/='CSF_SF'.and.dm_basis/='SF_SO'.and.&
            dm_basis/='CSF_SO'.and.dm_basis/='ALL') then
          call WarningMessage(2,'Unknown option for DMBAsis')
          call abend()
        endif
      case ('DIPO')
        flag_dipole=.True.
      case ('EMIS')
        flag_emiss=.True.
! reading pulse characteristic
      case ('NPUL')
        read(luin,'(I8)') N_pulse
        if (N_pulse>1) then
          call mma_deallocate(amp)
          call mma_deallocate(taushift)
          call mma_deallocate(pulse_vector)
          call mma_deallocate(sigma)
          call mma_deallocate(omega)
          call mma_deallocate(phi)
          call mma_allocate(taushift,N_pulse)
          call mma_allocate(amp,N_pulse)
          call mma_allocate(pulse_vector,N_pulse,3)
          call mma_allocate(sigma,N_pulse)
          call mma_allocate(omega,N_pulse)
          call mma_allocate(phi,N_pulse)
          do i=1,N_pulse
            amp(i)           = 0.0d0
            taushift(i)      = 0.0d0
            pulse_vector(i,1)= one
            pulse_vector(i,2)= zero
            pulse_vector(i,3)= zero
            sigma(i)         = 1.0d0/auToFs
            omega(i)         = 710d0/autoev
            phi(i)           = 0d0*pi
          enddo
        else if (N_pulse==0) then
          flag_pulse=.False.
        endif
      case ('PTYP')
        read(luin,'(A)') pulse_type
        call upCase(pulse_type)
        if (pulse_type(1:3)/='COS'.and.&
            pulse_type(1:3)/='SIN'.and.&
            pulse_type/='GAUSS'.and.&
            pulse_type/='MONO'.and.&
            pulse_type/='MONO_R_CIRCLE'.and.&
            pulse_type/='GAUSS_L_CIRCLE'.and.&
            pulse_type/='GAUSS_R_CIRCLE'.and.&
            pulse_type/='MONO_L_CIRCLE') then
          call WarningMessage(2,'Unknown option for PTYPe')
          call abend()
        endif
        ! check which power has the shape in case of SIN^N or COS^N
        if (pulse_type(1:3)=='COS'.or.&
            pulse_type(1:3)=='SIN') then
          read(pulse_type(4:),*) power_shape
        endif
      case ('TAUS')
        read(luin,*) (taushift(i),i=1,N_pulse)
        taushift=taushift/auToFs
      case ('AMPL')
        read(luin,*) (amp(i),i=1,N_pulse)
      case ('POLA')
        do i=1,N_pulse
          read(luin,*) (pulse_vector(i,j),j=1,3)
        enddo
      case ('SIGM')
        read(luin,*) (sigma(i),i=1,N_pulse)
        sigma=sigma/auToFs
        ! old value of sigma was underestimated by factor of pi/2
      case ('OMEG')
        read(luin,*) (omega(i),i=1,N_pulse)
        omega=omega/autoev
      case ('PHAS')
        read(luin,*) (phi(i),i=1,N_pulse)
        phi=phi*pi
      case ('CHIR')
        read(luin,*) linear_chirp
      case ('ACOR')
        flag_acorrection = .True.
      case('END ')
        exit
      case default
        write(u6,*) 'The corresponding keyword: ',line,' is unknown!'
        call abend()
    end select
  enddo
300   continue

  if (ipglob>1) then
    call dashes()
    write(u6,*)        'Input variables '
    call dashes()
    write(u6,sint)     'Number of spin manifolds:', N
    call dashes()
    write(u6,*) '      N       DET      CSF     STATES     SPIN'
    do i=1,N
      write(u6,'(5(i8,x))') i,ndet(i),nconf(i),lroots(i),ispin(i)
    enddo
    call dashes()
    write(u6,scha) 'State basis to be populated:', trim(p_style)
    write(u6,sint) 'Number of populated states:', n_populated
    if (p_style=='SO_THERMAL'.or.p_style=='SF_THERMAL') then
      write(u6,sdbl)   'Temperature:',    T
    endif
    write(u6,scha)     'Basis for propagation:',trim(basis)
    write(u6,sint)     'Number of states:',Nstate
    write(u6,sdbl)     'Initial time:',   initialtime*auToFs
    write(u6,sdbl)     'Final time:',     finaltime*auToFs
    write(u6,slog)     'Auger Decay:',    flag_decay
    write(u6,slog)     'Dissipation:',    flag_diss
    write(u6,slog)     'Ionization: ',    flag_dyson
    write(u6,slog)     'Pulse:',          flag_pulse
    if (flag_diss) then
!          write(6,sdbl)   'DeltaE:',         deltaE
!          write(6,sdbl)   'Coupling (cm-1):',V
      write(u6,sdbl)   'Gamma (Hartree):', gamma
    endif
    call dashes()
    write(u6,*)        'Pulse characteristics:'
    write(u6,sint) 'Number of pulses:',    N_pulse
    call dashes()
    if (flag_pulse.and.amp(1)/=0d0) then
      do i=1,N_pulse
        write(u6,scha) 'Pulse type:',     trim(pulse_type)
        write(u6,sint) 'Pulse # ',        i
        write(u6,sdbl) 'Amp:',            amp(i)
        write(u6,sdbl) 'Center:',         taushift(i)*auToFs
        write(u6,scmp) 'Polarization x:', pulse_vector(i,1)
        write(u6,scmp) 'Polarization y:', pulse_vector(i,2)
        write(u6,scmp) 'Polarization z:', pulse_vector(i,3)
      enddo
    endif
    call dashes()
  endif
  call close_luSpool(luin)
!  call qExit('read_input')! deprecated
  return
end
