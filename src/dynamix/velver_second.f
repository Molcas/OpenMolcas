************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2006, Igor Schapiro                                    *
************************************************************************
C
C *********************************************************************
C *                                                                   *
C * Second part of the velocity Verlet algorithm, which calculates    *
C * velocities of the next time step from the velocities from the     *
C * half step and the new forces.                                     *
C *                                                                   *
C * The algorithm is based on the example F4 of molecular dynamics    *
C * bible from Allen and Tildesley.                                   *
C *                                                                   *
C * 15/10/2006                                                        *
C * Igor Schapiro                                                     *
C *                                                                   *
C *********************************************************************

C REFERENCE:
C SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE VelVer_Second(irc)
      USE Isotopes
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
      INTEGER Iso
C
C     The parameter conv converts the gradients (Hartree/Bohr) to
C     (i)  forces (Hartree/Bohr)    => -1.0d0
C     (ii) forces (kJ/mole/Agstrom) => -4961.475514610d0
C
      PARAMETER   (conv=-1.0d0)
c     PARAMETER   (conv=-CONV_AU_TO_KJ_PER_MOLE_/Angstrom)
      PARAMETER   (kb = CONST_BOLTZMANN_/
     &             (CONV_AU_TO_KJ_*1.0D3))

*
      IF(IPRINT.EQ.INSANE) WRITE(6,*)' Entering ',ROUTINE
      CALL QENTER(ROUTINE)

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
      CALL Get_nAtoms_All(matom)
      CALL Get_Mass_All(Mass,matom)

C Check if reduced dimensionality
      IF (POUT .NE. 0) THEN
        CALL project_out(vel,force,natom)
      ENDIF
C
C     Definition of the time step
C
      DT2 = DT / 2.0d0
      Ekin = 0.0d0

      CALL Get_dScalar('MD_Time',time)

      DO i=1, natom
C     Determines the mass of an atom from its name
         IF (i.GT.matom) THEN
            CALL LeftAd(atom(i))
            Iso=0
            CALL Isotope(Iso,atom(i),Mass(i))
         END IF
C-------------------------------------------
         DO j=1, 3
            vel(3*(i-1)+j) = vel(3*(i-1)+j) + DT2 * force(3*(i-1)+j) /
     &      Mass(i)
         END DO
      END DO
C  Calling the thermostats for canonical ensemble
      IF (THERMO.eq.2) THEN
         CALL NhcThermo(vel)
      END IF
C

C Final kinetic energy
      DO i=1, natom
         DO j=1, 3
              Ekin = Ekin + 5.0d-01 * Mass(i) * (vel(3*(i-1)+j) ** 2)
         END DO
      END DO

      Call Add_Info('EKin',[EKin],1,6)
C

C
C     Write out the velocities
C
      caption = 'Velocities'
      lastline = ''
      CALL DxPtTableCo(caption,time,natom,atom,vel,lastline,Mass,force)
C
C     Read in the potential energy
C
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
C
C     Output
C
      IF (hybrid) THEN
          WRITE(6,'(//,5X,A,F8.1)') 'Final QM/MM Energy at time ',time
      ELSE
          WRITE(6,'(//,5X,A,F8.1)') 'Final Energy at time ',time
      ENDIF
      WRITE(6,'(5X,A,/)') '============================'
      WRITE(6,'(5X,A,6X,D19.12,1X,A)') 'Kinetic energy',Ekin,'a.u.'
      WRITE(6,'(5X,A,4X,D19.12,1X,A)') 'Potential Energy',Epot,'a.u.'
      WRITE(6,'(5X,A,8X,D19.12,1X,A)') 'Total Energy',Etot,'a.u.'

C--------------------------------------------------------------------C
C CANONICAL ENSEMBLE
C--------------------------------------------------------------------C
      tempNow = 2.D0 * Ekin /(3.D0*dble(natom)*kb)
      IF (THERMO.eq.2) THEN
         tempNow = 2.D0 * Ekin /(3.D0*dble(natom)*kb)
         WRITE(6,'(//,5X,A)') 'Canonical Ensemble'
         WRITE(6,'(5X,A)') 'The temperature is control with a '
         WRITE(6,'(5X,A)') 'Nose-Hoover chain of thermostats'
         WRITE(6,'(5X,A,/)') '========================'
         WRITE(6,'(5X,A,5X,E11.4,1X,A,//)') 'instantaneous temperature'
     &                                  ,tempNow,'kelvin'

      ENDIF
C--------------------------------------------------------------------C

C---------------------------------------------------------------------C
C     MICRO-CANONICAL ENSEMBLE
C---------------------------------------------------------------------C
      IF (THERMO.EQ.1) THEN
         CALL Get_dScalar('MD_Etot0',Etot0)
         WRITE(6,'(//,5X,A)') 'Micro-Canonical Ensemble'
         WRITE(6,'(5X,A,/)') '========================'
         WRITE(6,'(5X,A,7X,D11.4,1X,A)') 'Target Total Energy'
     &                                  ,Etot0,'a.u.'
         WRITE(6,'(5X,A,D11.4,1X,A)') 'Deviation from this Energy'
     &                              ,ABS(Etot0-Etot),'a.u.'
C
C     Check if the total energy is conserved and scale the velocities
C     if necessary.
C
C        1.0K * k_B
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
C---------------------------------------------------------------------C
      CALL Put_dScalar('MD_Time',time)

      Etot=Epot+Ekin

C------ Tully reascaling velocities in case of HOP ---------------------C
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
C         write (6,*) 'DEBUGS'
C         write (6,*) Epot, EtotLastPoint, Etot
C         write (6,*) 'EtotLastPoint'
C         write (6,*) EtotLastPoint
C         write (6,*) 'scalfac'
C         write (6,*) scalfac
C         write (6,*) 'nsAtom'
C         write (6,*) nsAtom
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

C
C     The return code is set in order to continue the loop
C
      irc=_RC_ALL_IS_WELL_
      CALL qExit(ROUTINE)
      RETURN
      END
