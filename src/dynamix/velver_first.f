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
C * First part of the velocity Verlet algorithm, which calculates the *
C * new positions for the next time step and the new velocities for   *
C * a half time step. The algorithm is based on the example F4 from   *
C * molecular dynamics bible of Allen and Tildesley.                  *
C *                                                                   *
C * 15/10/2006                                                        *
C * Igor Schapiro                                                     *
C *                                                                   *
C *********************************************************************

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE VelVer_First(irc)
      USE Isotopes
#include "prgm.fh"
#include "warnings.fh"
#include "Molcas.fh"
      PARAMETER   (ROUTINE='VV_First')
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
      EXTERNAL    IsFreeUnit
      INTEGER     natom,i,j,irc,file,IsFreeUnit,ipCoord
      REAL*8      DT2,DTSQ2,Ekin,time,totimpl,RMS

      CHARACTER  caption*15, lastline*80, filname*80
      LOGICAL    hybrid,qmmm
*
      INTEGER    natom2
*
      REAL*8, ALLOCATABLE ::      vel(:),xyz(:),force(:)
      REAL*8, ALLOCATABLE ::      Mass(:),tstxyz(:)
      CHARACTER, ALLOCATABLE ::  atom(:)*2, atom2(:)*2
      REAL*8, ALLOCATABLE ::     force2(:),xyz2(:)
      INTEGER Iso
*
      IF(IPRINT.EQ.INSANE) WRITE(6,*)' Entering ',ROUTINE
      CALL QENTER(ROUTINE)

      WRITE(6,*)'*** First step of the Velocity Verlet algorithm ***'
C
C     Check for QM/MM calculation
C
      filname = 'comqum.dat'
      CALL F_INQUIRE(filname,hybrid)
      hybrid=.False.
C
C     Read atom, their coordinates and forces
C
      IF (hybrid) THEN
         WRITE(6,'(/,5X,A)') 'Perform QM/MM Molecular Dynamics'
         CALL DxRdNAtomHbrd(natom)
      ELSE
         CALL DxRdNAtomStnd(natom)
      END IF

      CALL mma_allocate(vel,natom*3)
      CALL mma_allocate(xyz,natom*3)
      CALL mma_allocate(force,natom*3)
      CALL mma_allocate(tstxyz,natom*3)
      CALL mma_allocate(Mass,natom)
      CALL mma_allocate(atom,natom)

      IF (hybrid) THEN
         CALL DxRdHbrd(natom,atom,xyz,force)
      ELSE
         CALL DxRdStnd(natom,atom,xyz,force)
      END IF
      CALL Get_dScalar('MD_time',time)
C
C     Read the velocities
C
      CALL Get_Velocity(vel,3*natom)
C
C     Initialize the Mass variable
      CALL Get_nAtoms_All(matom)
      CALL Get_Mass_All(Mass,matom)
C
C Check if reduced dimensionality
      IF (POUT .NE. 0) THEN
        CALL project_out(vel,force,natom)
      ENDIF
C--------------------------------------------------------------------C
C CANONICAL ENSEMBLE
C--------------------------------------------------------------------C
      IF (THERMO.eq.2) THEN
         call NhcThermo(vel)
      ENDIF
C--------------------------------------------------------------------C

C     Write out the old coordinates
C
      CALL DxRdNAtomStnd(natom2)
      CALL mma_allocate(xyz2,natom*3)
      CALL mma_allocate(force2,natom*3)
      CALL mma_allocate(atom2,natom)
      CALL DxRdStnd(natom2,atom2,xyz2,force2)
**      IF (iPrint.ge.VERBOSE) THEN
         caption = 'Old Coordinates'
         lastline = ''
         CALL DxPtTableCo(caption,time,natom2,atom2,xyz2,lastline,
     &        Mass,force)
**      END IF
C
C     Definition of the time step
C
      DT2   = DT / 2.0D0
      DTSQ2 = DT * DT2
*
      Ekin = 0.0D0
      RMS = 0.0D0
      totimpl = 0.0D0
*
      DO i=1, natom
C     Determines the mass of an atom from its name
        IF (i.GT.matom) THEN
           CALL LeftAd(atom(i))
           Iso=0
           CALL Isotope(Iso,atom(i),Mass(i))
        END IF
C-------------------------------------------
        DO j=1, 3
*** Root mean square deviation ****************************
          tstxyz(3*(i-1)+j) = xyz(3*(i-1)+j)
***********************************************************
          xyz(3*(i-1)+j) = xyz(3*(i-1)+j) + DT *
     &    vel(3*(i-1)+j) + DTSQ2 * force(3*(i-1)+j) / Mass(i)
***********************************************************
          RMS=RMS+(tstxyz(3*(i-1)+j)-xyz(3*(i-1)+j))**2
***********************************************************
          Ekin = Ekin + 5.0D-01 * Mass(i) * (vel(3*(i-1)+j) ** 2)
          vel(3*(i-1)+j) = vel(3*(i-1)+j) + DT2 * force(3*(i-1)+j) /
     &       Mass(i)
          totimpl = totimpl + vel(3*(i-1)+j) * Mass(i)
        END DO
      END DO

      Call Add_Info('EKin',[EKin],1,6)

      RMS=SQRT(RMS/natom)
C
C     Update the Link-Atom position
C     (Temporary solution uses ipCoord instead of xyz)
C
      qmmm = .False.
         Call DecideOnESPF(qmmm)
         If (qmmm) Then
            Call GetMem('Coordinates','ALLO','REAL',ipCoord,3*natom)
            call dcopy_(3*natom,xyz,1,Work(ipCoord),1)
            Call LA_Morok(natom,ipCoord,2)
            call dcopy_(3*natom,Work(ipCoord),1,xyz,1)
            Call Free_Work(ipCoord)
         End If
C
C     Output
C
      IF (iPrint.ge.USUAL) THEN
        WRITE(6,400) 'Molecular Dynamics specifications (time = ',
     &                   time,' a.u.)'
        WRITE(6,'(5X,A,/)') '=========================================='
        WRITE(6,402) 'Kinetic energy',Ekin,'a.u.'
        WRITE(6,405) 'Total linear momentum ',totimpl,'a.u.'
        WRITE(6,402) 'RMS deviation ',RMS,'Bohr'
      END IF

      time = time + DT
C
C     Set the initialization flag and save the current time-value
C
      CALL Put_dScalar('MD_Time',time)
#ifdef _HDF5_
      call mh5_put_dset(dyn_time,time)
#endif
C
C     Write out the new coordinates
C
      caption = 'New Coordinates'
      lastline = ''
      CALL DxPtTableCo(caption,time,natom,atom,xyz,lastline,Mass,force)
C
C     Write coordinates to output file
C
*#ifdef _DEBUG_
*      WRITE(6,*)' Dynamix calls 2 DxCoord.'
*#endif
      CALL DxCoord(natom,atom,xyz,hybrid)
*#ifdef _DEBUG_
*      WRITE(6,*)' Dynamix back from 2 DxCoord.'
*#endif
C
C     Save the new coordinates
C
      IF (hybrid) THEN
          call dscal_(natom*3,Angstrom,xyz,1)
          file=IsFreeUnit(81)
          filname='prmcrd2'
          Call Molcas_Open(file,filname)
          WRITE(file,403)natom
          WRITE(file,404)(xyz(i),i=1,natom*3)
          CLOSE(file)
      ELSE
          CALL Put_Coord_Full(xyz,natom)
#ifdef _HDF5_
          call mh5_put_dset(dyn_geom,xyz)
#endif
      END IF

      CALL Put_Velocity(vel,3*natom)
#ifdef _HDF5_
      call mh5_put_dset(dyn_vel,vel)
#endif

      CALL mma_deallocate(vel)
      CALL mma_deallocate(xyz)
      CALL mma_deallocate(force)
      CALL mma_deallocate(tstxyz)
      CALL mma_deallocate(Mass)
      CALL mma_deallocate(atom)
      CALL mma_deallocate(xyz2)
      CALL mma_deallocate(force2)
      CALL mma_deallocate(atom2)
C
C     The return code is set in order to continue the loop
C
      irc=_RC_ALL_IS_WELL_
C
 400  FORMAT(5X,A,F8.1,A)
 402  FORMAT(5X,A14,8X,D11.4,1X,A)
 403  FORMAT(/,I5)
 404  FORMAT(6F12.7)
 405  FORMAT(5X,A22,D11.4,1X,A)
*
      CALL qExit(ROUTINE)
      RETURN
      END
