subroutine read_input()
  use rhodyn_data
  use rhodyn_utils, only: dashes
  implicit none
!
!
!
  integer              :: luin
  character(len=8)     :: input_id = '&RHODYN', tryname
  character(len=256)   :: line

  call qEnter('read_input')

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
        tout=tout*fstoau
      case('INIT')
        read(luin,*) initialtime
        initialtime=initialtime*fstoau
      case('FINA')
        read (luin,*) finaltime
        finaltime=finaltime*fstoau
      case('TSTE')
        read(luin,*) timestep
        timestep=timestep*fstoau
      case('METH')
        read(luin,*) method
      case('RK45')
        read(luin,*) errorthreshold
      case('RKSA')
        read(luin,*) safety
      case('DELT')
        read(luin,*) deltaE
        deltaE=deltaE*cmtoau
      case ('VCOU')
        read (luin,*) V
        V=V*cmtoau
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
!     tau_L3=tau_L3*fstoau
!     tau_L2=tau_L2*fstoau
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
        gamma=gamma*cmtoau
      case ('HRSO')
        HRSO=.True.
			case ('KEXT')
			  kext=.True.
      case ('IFPU')
        flag_pulse=.False.
      case ('TFDM')
        read(luin,*) time_fdm
        time_fdm=time_fdm*fstoau
        flag_fdm = .True.
      case ('DMBA')
        read(luin,'(A)') dm_basis
		  case ('DIPO')
        flag_dipole=.True.
      case ('EMIS')
        flag_emiss=.True.
      case ('PTYP')
        read(luin,'(A)') pulse_type
      case ('NPUL')
        read(luin,'(I8)') N_pulse
        if (N_pulse/=1) then
         deallocate(shift)
         deallocate(amp)
         deallocate(pulse_vector)
         deallocate(sigma)
         deallocate(omega)
         deallocate(phi)
         allocate(shift(N_pulse))
         allocate(amp(N_pulse))
         allocate(pulse_vector(N_pulse,3))
         allocate(sigma(N_pulse))
         allocate(omega(N_pulse))
         allocate(phi(N_pulse))
          do i=1,N_pulse
            amp(i)           = 2.5d0
            pulse_vector(i,1)= one
            pulse_vector(i,2)= zero
            pulse_vector(i,3)= zero
            sigma(i)         = autoev/5d0
            omega(i)         = 710d0/autoev
            Phi(i)           = 0d0*pi
          enddo
        endif
      case ('SHIF')
        shift = 0
        read(luin,*) (shift(i),i=2,N_pulse)
        shift=shift*fstoau
      case ('STST')
        read(luin,*) sin_tstar
      case ('STEN')
        read(luin,*) sin_tend
      case ('SSCA')
        read(luin,*) sin_scal
      case ('AMPL')
        read(luin,*) (amp(i),i=1,N_pulse)
      case ('TAUS')
        read(luin,*) tau
        tau =tau*fstoau
      case ('POLA')
        do i=1,N_pulse
          read(luin,*) (pulse_vector(i,j),j=1,3)
        enddo
      case ('SIGM')
        read(luin,*) (sigma(i),i=1,N_pulse)
        sigma=autoev/sigma
      case ('OMEG')
        read(luin,*) (omega(i),i=1,N_pulse)
        omega=omega/autoev
      case ('PHAS')
        read(luin,*) (phi(i),i=1,N_pulse)
        phi=phi*pi
      case('END ')
        exit
      case default
        print*,'The corresponding keyword: ',line,' is unknown!'
        stop
    end select
  enddo
300   continue

  if (ipglob>1) then
    call dashes()
    write(6,*)        'Input variables '
    call dashes()
    write(6,sint)     'Number of spin manifolds:', N
    call dashes()
    write(6,*) '      N       DET      CSF     STATES     SPIN'
    do i=1,N
      write(6,'(5(i8,x))') i,ndet(i),nconf(i),lroots(i),ispin(i)
    enddo
    call dashes()
    write(6,scha) 'State basis to be populated:', trim(p_style)
    write(6,sint) 'Number of populated states:', n_populated
    if (p_style=='SO_THERMAL'.or.p_style=='SF_THERMAL') then
      write(6,sdbl)   'Temperature:',    T
    endif
    write(6,scha)     'Basis for propagation:',trim(basis)
    write(6,sint)     'Number of states:',Nstate
    write(6,sdbl)     'Initial time:',   initialtime/fstoau
    write(6,sdbl)     'Final time:',     finaltime/fstoau
    write(6,slog)     'Auger Decay:',    flag_decay
    write(6,slog)     'Dissipation:',    flag_diss
    write(6,slog)     'Ionization: ',    flag_dyson
    write(6,slog)     'Pulse:',          flag_pulse
    if (flag_diss) then
!          write(6,sdbl)   'DeltaE:',         deltaE
!          write(6,sdbl)   'Coupling (cm-1):',V
      write(6,sdbl)   'Gamma (Hartree):', gamma
    endif
    call dashes()
    write(6,*)        'Pulse characteristics:'
    call dashes()
    if (flag_pulse.and.amp(1)/=0d0) then
      write(6,scha)   'Pulse type:',     trim(pulse_type)
      if (pulse_type=='sine_square') then
        write(6,sdbl) 'sin_tstar',       sin_tstar
        write(6,sdbl) 'sin_tend',        sin_tend
        write(6,sdbl) 'sin_scal',        sin_scal
      elseif (pulse_type=='Gaussian') then
        write(6,sdbl) 'Amp:',            amp(1)
        write(6,sdbl) 'Tau:',            tau/fstoau
        write(6,scmp) 'Polarization x:', pulse_vector(1,1)
        write(6,scmp) 'Polarization y:', pulse_vector(1,2)
        write(6,scmp) 'Polarization z:', pulse_vector(1,3)
      elseif (pulse_type=='Train_Diff'.or.pulse_type=='Train_Same') then
        write(6,sint) '# of pulses:',    N_pulse
        do i=1,N_pulse
          write(6,sint) 'Pulse # ',        i
          write(6,sdbl) 'Amp:',            amp(i)
          write(6,sdbl) 'Center:',      (tau+(i-1)*shift(i))/fstoau
          write(6,scmp) 'Polarization x:', pulse_vector(i,1)
          write(6,scmp) 'Polarization y:', pulse_vector(i,2)
          write(6,scmp) 'Polarization z:', pulse_vector(i,3)
        enddo
      endif
    endif
    call dashes()
  endif
  call close_luSpool(luin)
  call qExit('read_input')
  return
end
