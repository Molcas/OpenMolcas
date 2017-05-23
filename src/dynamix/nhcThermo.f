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
* Copyright (C) 2012, Felipe Zapata                                    *
************************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The Nose-Hoover chain of thermostat is based on the following paper
C Journal of Physical Chemistry B, 2001, 105, 7598
C The implemetation required first initialized the thermostat, then
C before calling the first part of velocity verlet the NHC is call,
C then the positions and velocities are updated using the verlet_first
C and after that the verlet_second is called and after this step
C the NHC subroutine is call. The order must not be changed
C
C Felipe Zapata, November 1, 2012
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
      SUBROUTINE NhcThermo (vel)
      USE Isotopes
      IMPLICIT REAL*8 (a-h,o-z)
#include "prgm.fh"
#include "Molcas.fh"
      PARAMETER    (ROUTINE='NhcThermo')
#include "warnings.fh"
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
      PARAMETER  (nh=6)
      INTEGER     natom,natom3,i,j,nIsoAtoms
      REAL*8      Ekin,kb
      PARAMETER   (kb = CONST_BOLTZMANN_/
     &             CONV_AU_TO_KJ_*1.0D3)
      REAL*8      NHC(nh), vel(*)
      LOGICAL     lIsotope
      CHARACTER, ALLOCATABLE ::   atom(:)*2
      REAL*8, ALLOCATABLE ::      Mass(:),dIsotopes(:)
      INTEGER Iso

*nh stands for the number of variables in the thermostat
* NHC = Q1,Q2,X1,X2,VX1,VX2,Scale
*
C
C    The parameter kb is the Boltzmann constant in a.u

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

CC READ PARAMETERS FROM RUNFILE

      CALL Get_nAtoms_Full(natom)

      CALL mma_allocate(atom,natom)
      CALL mma_allocate(Mass,natom)

      CALL Get_Name_Full(atom)

      natom3=3*natom

C     Read Thermostat Variables
      CALL Get_NHC(NHC,nh)

      Q1 = NHC(1)
      Q2 = NHC(2)
      X1 = NHC(3)
      X2 = NHC(4)
      Vx1 = NHC(5)
      Vx2 = NHC(6)

C     Check if the Isotope option is selected
C
      CALL Qpg_dArray('Isotopes',lIsotope,nIsoAtoms)
      IF (lIsotope) THEN
         CALL mma_allocate(dIsotopes,natom)
         CALL Get_dArray('Isotopes',dIsotopes,natom)
      WRITE(6,*) 'Isotopes Label:' ,lIsotope
      END IF

C     Initialize the Mass variable
C
      DO i=1, natom
         Mass(i)=0.D0
      END DO

      Ekin = 0.0D0


      DO i=1, natom

        CALL LeftAd(atom(i))
        Iso=0
        CALL Isotope(Iso,atom(i),Mass(i))
C     Manual isotope modification -----------
        IF (lIsotope) THEN
           IF ((dIsotopes(i)).NE.0.0D0) THEN
              Mass(i)=dIsotopes(i)
           END IF
        END IF
C-------------------------------------------
        DO j=1, 3
          Ekin = Ekin + 5.0D-01 * Mass(i) * (vel(3*(i-1)+j) ** 2.D0)
        END DO
      END DO



C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      DT2   = DT * 0.5D0
      DT4   = DT * 2.5D-1
      DT8   = DT * 1.25D-1

      G2 = (Q1*Vx1*Vx1 - TEMP*kb )/Q2
      Vx2 = Vx2 + G2*DT4
      Vx1 = Vx1 * exp(-Vx2*DT8)
      G1 = (2.D0*Ekin - 3.D0*dble(natom)*TEMP*kb)/Q1
      Vx1 = Vx1 + G1*DT4
      Vx1 = Vx1 *exp(-Vx2*DT8)
      sc  = exp(-Vx1*DT2)

      DO i=1, natom
          DO j=1, 3
             vel(3*(i-1)+j) = vel(3*(i-1)+j)*sc
          END DO
      END DO

      Ekin = Ekin*sc*sc

      X1 = X1 + Vx1*DT2
      X2 = X2 + Vx2*DT2
      Vx1 = Vx1 * exp(-Vx2*DT8)
      G1 = (2.D0*Ekin - 3.D0*dble(natom)*TEMP*kb)/Q1
      Vx1  = Vx1 + G1*DT4
      Vx1 = Vx1*exp(-Vx2*DT8)
      G2 = (Q1*(Vx1**2.D0)- TEMP*kb )/Q2
      Vx2 = Vx2 + G2*DT4

      NHC(3) = X1
      NHC(4) = X2
      NHC(5) = Vx1
      NHC(6) = Vx2

      CALL Put_NHC(NHC,nh)
#ifdef _HDF5_
      call mh5_put_dset(dyn_nh,NHC)
#endif

      CALL mma_deallocate(atom)
      CALL mma_deallocate(Mass)

      RETURN

      END
