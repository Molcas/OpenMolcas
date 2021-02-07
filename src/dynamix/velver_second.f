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
! Copyright (C) 2006, Igor Schapiro                                    *
!***********************************************************************
!
! *********************************************************************
! *                                                                   *
! * Second part of the velocity Verlet algorithm, which calculates    *
! * velocities of the next time step from the velocities from the     *
! * half step and the new forces.                                     *
! *                                                                   *
! * The algorithm is based on the example F4 of molecular dynamics    *
! * bible from Allen and Tildesley.                                   *
! *                                                                   *
! * 15/10/2006                                                        *
! * Igor Schapiro                                                     *
! *                                                                   *
! *********************************************************************

! REFERENCE:
! SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.

!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE VelVer_Second(irc)
#ifdef _HDF5_
      USE mh5, ONLY: mh5_put_dset
#endif
      IMPLICIT REAL*8 (a-h,o-z)
#include "prgm.fh"
#include "warnings.fh"
#include "Molcas.fh"
      PARAMETER    (ROUTINE='VV_Second')
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
      EXTERNAL     IsFreeUnit
      INTEGER      i,j,natom,file,IsFreeUnit,nsAtom
      CHARACTER    filname*80,line*80
      REAL*8       conv,tolerance,DT2,time,kb
      REAL*8       Ekin,Epot,Etot,Etot0,Ekin_target
      REAL*8       EtotLastPoint
      LOGICAL      hybrid,found
      CHARACTER  caption*15, lastline*80
      CHARACTER, ALLOCATABLE ::    atom(:)*2
      REAL*8, ALLOCATABLE ::       Mass(:),vel(:),force(:),xyz(:)
!
!     The parameter conv converts the gradients (Hartree/Bohr) to
!     (i)  forces (Hartree/Bohr)    => -1.0d0
!     (ii) forces (kJ/mole/Agstrom) => -4961.475514610d0
!
      PARAMETER   (conv=-1.0d0)
!     PARAMETER   (conv=-CONV_AU_TO_KJ_PER_MOLE_/Angstrom)
      PARAMETER   (kb = CONST_BOLTZMANN_/                               &
     &             (CONV_AU_TO_KJ_*1.0D3))

!
      IF(IPRINT.EQ.INSANE) WRITE(6,*)' Entering ',ROUTINE

      WRITE(6,*)'*** Second step of the Velocity Verlet algorithm ***'

      filname = 'comqum.dat'
      CALL F_INQUIRE(filname,hybrid)

      IF (hybrid) THEN
         WRITE(6,'(/,5X,A)') 'Perform QM/MM Molecular Dynamics'
         CALL DxRdNAtomHbrd(natom)
      ELSE
         CALL DxRdNAtomStnd(natom)
      END IF

      CALL mma_allocate(atom,natom)
      CALL mma_allocate(Mass,natom)
      CALL mma_allocate(vel,natom*3)
      CALL mma_allocate(force,natom*3)
      CALL mma_allocate(xyz,natom*3)

      IF (hybrid) THEN
         CALL DxRdHbrd(natom,atom,xyz,force)
      ELSE
         CALL DxRdStnd(natom,atom,xyz,force)
      END IF

      CALL Get_Velocity(vel,3*natom)
      CALL GetMassDx(Mass,natom)

! Check if reduced dimensionality
      IF (POUT .NE. 0) THEN
        CALL project_out_for(force,natom)
      ELSEIF (PIN .NE. natom*3) THEN
        CALL project_in_for(force,natom)
      ENDIF
!
!     Definition of the time step
!
      DT2 = DT / 2.0d0
      Ekin = 0.0d0

      CALL Get_dScalar('MD_Time',time)

      DO i=1, natom
         DO j=1, 3
            vel(3*(i-1)+j) = vel(3*(i-1)+j) + DT2 * force(3*(i-1)+j) /  &
     &      Mass(i)
         END DO
      END DO
!  Calling the thermostats for canonical ensemble
      IF (THERMO.eq.2) THEN
         CALL NhcThermo(vel)
      END IF
!
! Check if reduced dimensionality (should not be needed)
      IF (POUT .NE. 0) THEN
        CALL project_out_vel(vel,natom)
      ENDIF

! Final kinetic energy
      DO i=1, natom
         DO j=1, 3
              Ekin = Ekin + 5.0d-01 * Mass(i) * (vel(3*(i-1)+j) ** 2)
         END DO
      END DO

      Call Add_Info('EKin',[EKin],1,6)
!

!
!     Write out the velocities
!
      caption = 'Velocities'
      lastline = ''
      CALL DxPtTableCo(caption,time,natom,atom,vel,lastline,Mass,force)
!
!     Read in the potential energy
!
      IF (hybrid) THEN
          file=IsFreeUnit(81)
          filname = 'fixenergy.out'
          Call Molcas_Open(file,filname)
          READ(file,*)line
             DO WHILE (i.le.80)
                IF (line(i:i).eq.'$') THEN
                   READ(file,'(E20.13)')Epot
                   i=80
                ELSEIF (i.eq.80.and.line(i:i).ne.'$') THEN
                   WRITE(6,*)'No energy found'
                   CALL Abend()
                END IF
                i=i+1
             END DO
          CLOSE(file)
      ELSE
          CALL Get_dScalar('Last Energy',Epot)
      ENDIF

      Etot = Epot + Ekin
!
!     Output
!
      IF (hybrid) THEN
          WRITE(6,'(//,5X,A,F8.1)') 'Final QM/MM Energy at time ',time
      ELSE
          WRITE(6,'(//,5X,A,F8.1)') 'Final Energy at time ',time
      ENDIF
      WRITE(6,'(5X,A,/)') '============================'
      WRITE(6,'(5X,A,6X,D19.12,1X,A)') 'Kinetic energy',Ekin,'a.u.'
      WRITE(6,'(5X,A,4X,D19.12,1X,A)') 'Potential Energy',Epot,'a.u.'
      WRITE(6,'(5X,A,8X,D19.12,1X,A)') 'Total Energy',Etot,'a.u.'

!--------------------------------------------------------------------C
! CANONICAL ENSEMBLE
!--------------------------------------------------------------------C
      tempNow = 2.D0 * Ekin /(3.D0*dble(natom)*kb)
      IF (THERMO.eq.2) THEN
         tempNow = 2.D0 * Ekin /(3.D0*dble(natom)*kb)
         WRITE(6,'(//,5X,A)') 'Canonical Ensemble'
         WRITE(6,'(5X,A)') 'The temperature is control with a '
         WRITE(6,'(5X,A)') 'Nose-Hoover chain of thermostats'
         WRITE(6,'(5X,A,/)') '========================'
         WRITE(6,'(5X,A,5X,E11.4,1X,A,//)') 'instantaneous temperature' &
     &                                  ,tempNow,'kelvin'

      ENDIF
!--------------------------------------------------------------------C

!---------------------------------------------------------------------C
!     MICRO-CANONICAL ENSEMBLE
!---------------------------------------------------------------------C
      IF (THERMO.EQ.1) THEN
         CALL Get_dScalar('MD_Etot0',Etot0)
         WRITE(6,'(//,5X,A)') 'Micro-Canonical Ensemble'
         WRITE(6,'(5X,A,/)') '========================'
         WRITE(6,'(5X,A,7X,D11.4,1X,A)') 'Target Total Energy'          &
     &                                  ,Etot0,'a.u.'
         WRITE(6,'(5X,A,D11.4,1X,A)') 'Deviation from this Energy'      &
     &                              ,ABS(Etot0-Etot),'a.u.'
!
!     Check if the total energy is conserved and scale the velocities
!     if necessary.
!
!        1.0K * k_B
         tolerance=1.0D0*kb
         tolerance=1.5D0*natom*tolerance
         IF (ABS(Etot0-Etot).gt.tolerance) THEN
            Ekin_target=Etot0-Epot
            scalfac = sqrt(Ekin_target / Ekin)
            DO i=1, natom
               DO j=1,3
                  vel(3*(i-1)+j)=scalfac*vel(3*(i-1)+j)
               END DO
            END DO
            WRITE(6,'(5X,A)') 'is larger then Scaling-Threshold XX.'
            WRITE(6,'(5X,A)') 'Velocity scaling is necessary.'
            WRITE(6,401) 'Velocity scaling factor',scalfac
            Ekin=Ekin_target
         ELSE
            WRITE(6,'(5X,A)') 'is smaller then Scaling-Threshold XX.'
            WRITE(6,'(5X,A)') 'Velocity scaling is not necessary.'
         ENDIF
      ENDIF
!---------------------------------------------------------------------C
      CALL Put_dScalar('MD_Time',time)

      Etot=Epot+Ekin

!------ Tully reascaling velocities in case of HOP ---------------------C
      call qpg_iscalar('hopped',Found)
      if (found) call get_lscalar('hopped',Found)
      if (found) then
         call Get_iScalar('Unique atoms',nsAtom)
         call get_dScalar('MD_Etot',EtotLastPoint)
         Ekin_target=EtotLastPoint-Epot
         if (Ekin_target.lt.0.0) then
            Ekin_target=0.0
            write (6,*) 'warning negative kin energy rescaled to 0.0'
         end if
!         write (6,*) 'DEBUGS'
!         write (6,*) Epot, EtotLastPoint, Etot
!         write (6,*) 'EtotLastPoint'
!         write (6,*) EtotLastPoint
!         write (6,*) 'scalfac'
!         write (6,*) scalfac
!         write (6,*) 'nsAtom'
!         write (6,*) nsAtom
         scalfac = sqrt(Ekin_target / Ekin)
         write(6,*) 'Velocities before Hop:'
         do i=1, nsAtom
           write(6,*) vel(i*3-2), vel(i*3-1), vel(i*3)
         end do
         DO i=1, natom
            DO j=1,3
               vel(3*(i-1)+j)=scalfac*vel(3*(i-1)+j)
            END DO
         END DO

         write(6,*) 'Velocities after Hop:'
         do i=1, nsAtom
           write(6,*) vel(i*3-2), vel(i*3-1), vel(i*3)
         end do
         Etot=Epot+Ekin_target
      end if

      call Put_dScalar('MD_Etot',Etot)
#ifdef _HDF5_
       call mh5_put_dset(dyn_etot,Etot)
#endif
      CALL DxEnergies(time,Epot,Ekin,Etot)

      CALL DxWtVel(vel,3*natom)
      CALL Put_Velocity(vel,3*natom)
#ifdef _HDF5_
      call mh5_put_dset(dyn_vel,vel)
#endif

 401  FORMAT(5X,A,3X,E11.4)

      CALL mma_deallocate(atom)
      CALL mma_deallocate(Mass)
      CALL mma_deallocate(vel)
      CALL mma_deallocate(force)
      CALL mma_deallocate(xyz)

!
!     The return code is set in order to continue the loop
!
      irc=_RC_ALL_IS_WELL_
      RETURN
      END
