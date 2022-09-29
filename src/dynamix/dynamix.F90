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
use Dynamix_Globals, only: dyn_etot, dyn_etot0, dyn_fileid, dyn_nh, dyn_vel, File_H5Res, lH5Restart
#endif
use Dynamix_Globals, only: DT, PIN, POUT, THERMO, TEMP, RESTART, VELO, nh, iQ1, iQ2, iX1, iX2, iVx1, iVx2, VelVer, VV_First, &
                           VV_Second, Gromacs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: auTokJ, auTofs, kBoltzmann, Zero, One, Three, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp), parameter :: nTasks = 3
integer(kind=iwp) :: i, irc, iRlxRoot, iseed, iTask, Itr, j, LuInput, mTasks, MxItr, natom, nFlag, nRoots, Task(nTasks)
real(kind=wp) :: arg, buffer, Ekin, Epot, Etot0, Freq, mean, NHC(nh), Q1, Q2, Sigma, time, val
logical(kind=iwp) :: Found, lHop
character(len=16) :: StdIn
character(len=15) :: caption
character(len=8) :: ENV
real(kind=wp), parameter :: kb = kBoltzmann/(auTokJ*1.0e3_wp)
character(len=2), allocatable :: atom(:)
real(kind=wp), allocatable :: Mass(:), vel(:), pcoo(:,:)
integer(kind=iwp), external :: AixRm, IsFreeUnit, IsStructure
#include "warnings.h"

iReturn = 99

! Initialize Dynamix and set default values

#ifdef _DEBUGPRINT_
write(u6,*) ' Dynamix calls Init_Dynamix.'
#endif
call Init_Dynamix()
#ifdef _DEBUGPRINT_
write(u6,*) ' Dynamix back from Init_Dynamix.'
#endif

! Read the input

#ifdef _HDF5_
call cre_dyn()
#endif
#ifdef _DEBUGPRINT_
write(u6,*) ' Dynamix calls Readin_Dynamix.'
#endif
call Readin_Dynamix(Task,nTasks,mTasks)
#ifdef _DEBUGPRINT_
write(u6,*) ' Dynamix back from Readin_Dynamix.'
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
  if (RESTART == Zero) then
    time = Zero
  else
    time = RESTART
    write(u6,'(5X,A,T55,F9.2,A)') 'MD restart time = ',RESTART,' a.u.'
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
    Freq = One/(22.0_wp/auTofs)
    Q1 = Three*natom*TEMP*Kb/(Freq*Freq)
    Q2 = TEMP*Kb/(Freq*Freq)

    NHC(iQ1) = Q1
    NHC(iQ2) = Q2
    NHC(iX1) = Zero
    NHC(iX2) = Zero
    NHC(iVx1) = Zero
    NHC(iVx2) = Zero

    call Put_NHC(NHC,nh)
#   ifdef _HDF5_
    call mh5_put_dset(dyn_nh,NHC)
#   endif

  end if

  ! Check if nuclear coordinates to project out from the dynamics
  if ((POUT == 0) .and. (PIN == natom*3)) then
    write(u6,'(5X,A,T55)') 'Dynamics in full dimensionality.'
  else
    write(u6,'(5X,A,T55)') 'Dynamics in reduced dimensionality.'
    if (POUT /= 0) then
      call mma_allocate(pcoo,POUT,natom*3)
      call DxRdOut(pcoo,POUT,natom)
      ! Save on RUNFILE
      call Put_dArray('Proj_Coord',pcoo,POUT*natom*3)
    else if (PIN /= natom*3) then
      call mma_allocate(pcoo,PIN,natom*3)
      call DxRdIn(pcoo,PIN,natom)
      ! Save on RUNFILE
      call Put_dArray('Keep_Coord',pcoo,PIN*natom*3)
    end if
  end if

  if (VELO == 1) then
    call DxRdVel(vel,natom)
    write(u6,'(5X,A,T55)') 'The initial velocities (bohr/au) are read in.'
  else if (VELO == 2) then
    call DxRdVel(vel,natom)
    do i=1,natom
      do j=1,3
        vel(3*(i-1)+j) = vel(3*(i-1)+j)/sqrt(Mass(i))
      end do
    end do
    write(u6,'(5X,A,T55)') 'The initial mass weighted velocities (bohr/au) are read in.'

    ! Maxwell-Boltzmann distribution
  else if (VELO == 3) then
    nFlag = 0
    val = Zero
    buffer = Zero

    write(u6,'(5X,A,T55)') 'The initial velocities (bohr/au) are taken '
    write(u6,'(5X,A,f9.2,A)') 'from a Boltzmann distribution at',TEMP,' kelvin'

    call getSeed(iseed)

    do i=1,natom
      arg = TEMP*Kb/Mass(i)
      Sigma = sqrt(arg)
      mean = Zero
      do j=1,3
        call RandomGauss(mean,Sigma,iseed,nflag,buffer,Val)
        vel(3*(i-1)+j) = Val

        !write(u6,'(5x,a,t55,d16.8)') 'Vel = ',Val

      end do
    end do
  else
    do i=1,3*natom
      vel(i) = Zero
    end do
    write(u6,'(5X,A,T55)') 'The initial velocities are set to zero.'
  end if
  caption = 'Velocities'
  call DxPtTableWithoutMassForce(caption,time,natom,atom,vel)

  ! Check if reduced dimensionality
  if (POUT /= 0) then
    call project_out_vel(vel,natom)
  else if (PIN /= natom*3) then
    call project_in_vel(vel,natom)
    caption = 'Vel (red dim)'
    call DxPtTableWithoutMassForce(caption,time,natom,atom,vel)
  end if

  ! Calculate the kinetic energy
  if (VELO > 0) then
    Ekin = Zero
    do i=1,natom
      do j=1,3
        Ekin = Ekin+Half*Mass(i)*(vel(3*(i-1)+j)**2)
      end do
    end do
  else
    Ekin = Zero
  end if
  write(u6,'(5x,a,6x,d19.12,1x,a)') 'Kinetic energy',Ekin,'a.u.'
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
  write(u6,'(5x,a,8x,d19.12,1x,a)') 'Total Energy',Etot0,'a.u.'
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
      write(u6,*) ' Dynamix calls VelVer_Second.'
#     endif
      call VelVer_Second(irc)
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Dynamix back from VelVer_Second.'
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
        write(u6,*) ' Dynamix calls VelVer_First.'
#       endif
        call VelVer_First(irc)
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) ' Dynamix back from VelVer_First.'
#     endif
    else
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Dynamix calls VelVer_First.'
#     endif

      call VelVer_First(irc)
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Dynamix back from VelVer_First.'
#     endif
    end if
  else if (Task(iTask) == VV_First) then
    call VelVer_First(irc)

  else if (Task(iTask) == VV_Second) then
    call VelVer_Second(irc)

  else if (Task(iTask) == Gromacs) then
    call GROM(irc)

  else
    write(u6,*) 'Illegal task'
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
