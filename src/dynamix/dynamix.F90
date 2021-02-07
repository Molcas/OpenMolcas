!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Dynamix(iReturn)

#ifdef _HDF5_
use mh5, only: mh5_put_dset, mh5_close_file
#endif

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
integer AixRm
external IsFreeUnit, AixRm, IsStructure
parameter(nTasks=3)
parameter(nh=6)
character StdIn*16, caption*15, ENV*8
real*8 time, mean, kb
real*8 Epot, Ekin, Etot0
real*8 NHC(nh)
integer Task(nTasks), natom, IsFreeUnit, irc, Itr, MxItr
integer iRlxRoot, nRoots, i
logical Found, lHop
integer VelVer, VV_First, VV_Second, Gromacs, VV_Dump
parameter(au_time=CONST_AU_TIME_IN_SI_*1.0d15)
parameter(kb=CONST_BOLTZMANN_/(CONV_AU_TO_KJ_*1.0d3))
parameter(VelVer=1,VV_First=2,VV_Second=3,Gromacs=4,VV_Dump=5)
parameter(iQ1=1,iQ2=2,iX1=3,iX2=4,iVx1=5,iVx2=6)
character, allocatable :: atom(:)*2
real*8, allocatable :: Mass(:), vel(:), pcoo(:,:)

iReturn = 99

! Initialize Dynamix and set default values

#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix calls Init_Dynamix.'
#endif
call Init_Dynamix
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix back from Init_Dynamix.'
#endif

!     Read the input

#ifdef _HDF5_
call cre_dyn
#endif
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix calls Readin_Dynamix.'
#endif
call Readin_Dynamix(Task,nTasks,mTasks)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix back from Readin_Dynamix.'
#endif

! Check if this is an initial run of Dynamix

call Qpg_dScalar('MD_Time',Found)

#ifdef _HDF5_
if (.not. found .and. lH5Restart) then
  call restart_dynamix(file_h5res)
  found = .true.
end if
#endif

! Generate or read velocities if this is an initial run

if (.not. Found) then
  ! Check if the RESTART keyword was used.
  if (RESTART == 0.0d0) then
    time = 0.000d0
  else
    time = RESTART
    write(6,'(5X,A,T55,F9.2,A)') 'MD restart time = ',RESTART,' a.u.'
  end if

  call Put_dScalar('MD_Time',time)
  call Get_nAtoms_Full(natom)

  call mma_allocate(atom,natom)
  call mma_allocate(Mass,natom)
  call mma_allocate(vel,natom*3)

  call Get_Name_Full(atom)
  call GetMassDx(Mass,natom)

  ! Initialize Thermostat Variables

  if (THERMO == 2) then
    Freq = 1.d0/(2.2d1/au_time)
    Q1 = 3.d0*dble(natom)*TEMP*Kb/(Freq*Freq)
    Q2 = TEMP*Kb/(Freq*Freq)

    NHC(iQ1) = Q1
    NHC(iQ2) = Q2
    NHC(iX1) = 0.d0
    NHC(iX2) = 0.d0
    NHC(iVx1) = 0.d0
    NHC(iVx2) = 0.d0

    call Put_NHC(NHC,nh)
#   ifdef _HDF5_
    call mh5_put_dset(dyn_nh,NHC)
#   endif

  end if

  ! Check if nuclear coordinates to project out from the dynamics
  if ((POUT == 0) .and. (PIN == natom*3)) then
    write(6,'(5X,A,T55)') 'Dynamics in full dimensionality.'
  else
    write(6,'(5X,A,T55)') 'Dynamics in reduced dimensionality.'
    if (POUT /= 0) then
      call mma_allocate(pcoo,POUT,natom*3)
      call DxRdOut(pcoo,POUT,natom)
      ! Save on RUNFILE
      call Put_dArray('Proj_Coord',pcoo,POUT*natom*3)
    elseif (PIN /= natom*3) then
      call mma_allocate(pcoo,PIN,natom*3)
      call DxRdIn(pcoo,PIN,natom)
      ! Save on RUNFILE
      call Put_dArray('Keep_Coord',pcoo,PIN*natom*3)
    end if
  end if

  if (VELO == 1) then
    call DxRdVel(vel,natom)
    write(6,'(5X,A,T55)') 'The initial velocities (bohr/au) are read in.'
  elseif (VELO == 2) then
    call DxRdVel(vel,natom)
    do i=1,natom
      do j=1,3
        vel(3*(i-1)+j) = vel(3*(i-1)+j)/sqrt(Mass(i))
      end do
    end do
    write(6,'(5X,A,T55)') 'The initial mass weighted velocities (bohr/au) are read in.'

    ! Maxwell-Boltzmann distribution
  elseif (VELO == 3) then
    nFlag = 0
    val = 0.d0
    buffer = 0.d0

    write(6,'(5X,A,T55)') 'The initial velocities (bohr/au) are taken '
    write(6,'(5X,A,f9.2,A)') 'from a Boltzmann distribution at',TEMP,' kelvin'

    call getSeed(iseed)

    do i=1,natom
      arg = TEMP*Kb/Mass(i)
      Sigma = sqrt(arg)
      mean = 0.d0
      do j=1,3
        call RandomGauss(mean,Sigma,iseed,nflag,buffer,Val)
        vel(3*(i-1)+j) = Val

        !write(6,'(5X,A,T55,D16.8)') 'Vel = ',Val

      end do
    end do
  else
    do i=1,3*natom
      vel(i) = 0.000000000000d0
    end do
    write(6,'(5X,A,T55)') 'The initial velocities are set to zero.'
  end if
  caption = 'Velocities'
  call DxPtTableWithoutMassForce(caption,time,natom,atom,vel)

  ! Check if reduced dimensionality
  if (POUT /= 0) then
    call project_out_vel(vel,natom)
  elseif (PIN /= natom*3) then
    call project_in_vel(vel,natom)
    caption = 'Vel (red dim)'
    call DxPtTableWithoutMassForce(caption,time,natom,atom,vel)
  end if

  ! Calculate the kinetic energy
  if (VELO > 0) then
    Ekin = 0.000000000000d0
    do i=1,natom
      do j=1,3
        Ekin = Ekin+(5.0D-01)*Mass(i)*(vel(3*(i-1)+j)**2)
      end do
    end do
  else
    Ekin = 0.000000000000d0
  end if
  write(6,'(5X,A,6X,D19.12,1X,A)') 'Kinetic energy',Ekin,'a.u.'
  ! Save the velocities on RUNFILE
  call Put_Velocity(vel,3*natom)
  ! Save the total energy on RUNFILE if the total energy should be conserved.
  call Get_dScalar('Last Energy',Epot)
  Etot0 = Epot+Ekin
  call Put_dScalar('MD_Etot0',Etot0)
  call Put_dScalar('MD_Etot',Etot0)
# ifdef _HDF5_
  call mh5_put_dset(dyn_vel,vel)
  call mh5_put_dset(dyn_etot0,Etot0)
  call mh5_put_dset(dyn_etot,Etot0)
# endif
  call DxEnergies(time,Epot,Ekin,Etot0)
  write(6,'(5X,A,8X,D19.12,1X,A)') 'Total Energy',Etot0,'a.u.'
  call mma_deallocate(atom)
  call mma_deallocate(Mass)
  call mma_deallocate(vel)
  if ((POUT /= 0) .or. (PIN /= natom*3)) then
    call mma_deallocate(pcoo)
  end if
end if

! Execute the tasks

do iTask=1,mTasks

  if (Task(iTask) == VelVer) then

    if (Found) then

#     ifdef _DEBUGPRINT_
      write(6,*) ' Dynamix calls VelVer_Second.'
#     endif
      call VelVer_Second(irc)
#     ifdef _DEBUGPRINT_
      write(6,*) ' Dynamix back from VelVer_Second.'
#     endif

      ! Check for Hopping?

      lHop = .false.
      call qpg_iScalar('MaxHops',lHop)
      if (lHop) then

        ! Read the roots

        call Get_iScalar('Number of roots',nRoots)
        call Get_iScalar('Relax CASSCF root',iRlxRoot)

        ! Run RASSI

        LuInput = 11
        LuInput = IsFreeUnit(LuInput)
        call StdIn_Name(StdIn)
        call Molcas_Open(LuInput,StdIn)
        write(LuInput,'(A)') '>export DYN_OLD_TRAP=$MOLCAS_TRAP'
        write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'
        write(LuInput,'(A)') ' &RASSI &End'
        write(LuInput,'(A)') ' NR OF JOBIPHS'
        write(LuInput,*) ' 1 ',nRoots
        write(LuInput,*) (i,i=1,nRoots)
        !write(LuInput,'(X,I1,1X,I1)') inxtState,iRlxRoot
        write(LuInput,'(A)') ' HOP'
        write(LuInput,'(A)') 'End of Input'
        write(LuInput,'(A)') ' &Dynamix &End'
        write(LuInput,'(A)') ' VV_First'
        write(LuInput,'(A)') ' DT'
        write(LuInput,*) DT
        write(LuInput,'(A)') 'THERMO'
        write(LuInput,*) THERMO
        write(LuInput,'(A)') 'VELO'
        write(LuInput,*) VELO
        write(LuInput,'(A)') 'OUT'
        write(LuInput,*) POUT
        write(LuInput,'(A)') 'IN'
        write(LuInput,*) PIN
        write(LuInput,'(A)') 'End of Input'
        write(LuInput,'(A)') '>export MOLCAS_TRAP=$DYN_OLD_TRAP'
        close(LuInput)
        call Finish(_RC_INVOKED_OTHER_MODULE_)
      else
#       ifdef _DEBUGPRINT_
        write(6,*) ' Dynamix calls VelVer_First.'
#       endif
        call VelVer_First(irc)
      end if

#     ifdef _DEBUGPRINT_
      write(6,*) ' Dynamix back from VelVer_First.'
#     endif
    else
#     ifdef _DEBUGPRINT_
      write(6,*) ' Dynamix calls VelVer_First.'
#     endif

      call VelVer_First(irc)
#     ifdef _DEBUGPRINT_
      write(6,*) ' Dynamix back from VelVer_First.'
#     endif
    end if
  else if (Task(iTask) == VV_First) then
    call VelVer_First(irc)

  else if (Task(iTask) == VV_Second) then
    call VelVer_Second(irc)

  else if (Task(iTask) == Gromacs) then
    call GROM(irc)

  else if (Task(iTask) == VV_Dump) then
    call VelVer_Dump(irc)

  else
    write(6,*) 'Illegal task'
    call Abend()
  end if
end do

! Remove the GRADS file

call f_Inquire('GRADS',Found)
if (Found) then
  if (AixRm('GRADS') /= 0) call Abend()
end if

#ifdef _HDF5_
call mh5_close_file(dyn_fileid)
#endif

! If running in a DoWhile loop, we turn a successful
! return code into "continue loop", except on the
! last iteration

if ((IsStructure() == 1) .and. (irc == 0)) then
  MxItr = 0
  call GetEnvf('MOLCAS_MAXITER',ENV)
  if (ENV /= ' ') then
    read(ENV,*) MxItr
  end if
  Itr = 1
  call GetEnvf('MOLCAS_ITER',ENV)
  if (ENV /= ' ') then
    read(ENV,*) Itr
  end if
  if (Itr < MxItr) then
    iReturn = _RC_CONTINUE_LOOP_
  else
    iReturn = irc
  end if
else
  iReturn = irc
end if

return

end subroutine Dynamix
