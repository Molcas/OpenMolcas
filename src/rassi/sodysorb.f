************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE SODYSORB(NSS,LUTOTR,LUTOTI,DYSAMPS,SFDYS,NZ,
     &                 SODYSAMPS,SODYSAMPSR,SODYSAMPSI,SOENE)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "Files.fh"
      INTEGER   SOTOT,SFTOT,SO2SFNUM
      INTEGER   NZ,LSZZ,ORBNUM
      INTEGER   SODYSCIND
      INTEGER   INDJ,INDI,SFI,SFJ,ZI,ZJ,NSZZ,NDUM

      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION SFDYS(NZ,NSTATE,NSTATE)
      DIMENSION SODYSAMPS(NSS,NSS)
      DIMENSION SODYSAMPSR(NSS,NSS)
      DIMENSION SODYSAMPSI(NSS,NSS)
      DIMENSION SODYSCOFSR(NZ),SODYSCOFSI(NZ)
      DIMENSION SZZFULL(NZ,NZ)

      ! Joel
      DIMENSION SODYSCMOR(NZ*NZ)
      DIMENSION SODYSCMOI(NZ*NZ)
      DIMENSION SOENE(NZ)
      DIMENSION DYSEN(NZ)
      DIMENSION AMPS(NZ)
      INTEGER   IFILE

! +++ J.Norell 2018

C Calculates spin-orbit Dyson orbitals
C The routine was in some part adapted from DO_SONATORB

C 1. (Fast): Compute the SO dysamps from the SF dysamps and the
C      SO eigenvectors (approximation)
C 2. (Slower): Export Dyson orbitals in Molden and .DysOrb format
C      for all the requested initial states
C      The corresponding exact amplituds will be calculated.

! ****************************************************************

C Setup SO2SFNUM list which contains the original SF state numbers
C as a function of the SO state number

      SO2SFNUM=0
      CALL GETMEM('SO2SF','ALLO','INTE',SO2SFNUM,NSS)
      ISS=0
      SOTOT=0
      SFTOT=0
      DO ISTATE=1,NSTATE
       JOB1=iWork(lJBNUM+ISTATE-1)
       MPLET1=MLTPLT(JOB1)
       SFTOT=SFTOT+1

       DO MSPROJ1=-MPLET1+1,MPLET1-1,2
        SOTOT=SOTOT+1
        IWORK(SO2SFNUM+SOTOT-1)=SFTOT

       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE

      IF(.NOT.DYSEXPORT) THEN ! Approximative amplitude calculation

C Now construct the SF dysamps in the multiplicity expanded basis
C (initially all real, therefore put into SODYSAMPSR)
      SODYSAMPSR=0.0D0
      SODYSAMPSI=0.0D0
      DO JSTATE=1,NSS
         DO ISTATE=JSTATE+1,NSS
          SFJ=IWORK(SO2SFNUM+JSTATE-1)
          SFI=IWORK(SO2SFNUM+ISTATE-1)
          SODYSAMPSR(JSTATE,ISTATE)=DYSAMPS(SFJ,SFI)
          SODYSAMPSR(ISTATE,JSTATE)=DYSAMPS(SFJ,SFI)
         END DO
      END DO

C Now perform the transformation from SF dysamps to SO dysamps
C by combining the multiplicity expanded SF dysamps with the
C SO eigenvector in the ZTRNSF routine.
      CALL ZTRNSF(NSS,WORK(LUTOTR),WORK(LUTOTI),SODYSAMPSR,SODYSAMPSI)

C Compute the magnitude of the complex amplitudes
      SODYSAMPSR=SODYSAMPSR*SODYSAMPSR
      SODYSAMPSI=SODYSAMPSI*SODYSAMPSI
      SODYSAMPS=SQRT(SODYSAMPSR+SODYSAMPSI)

      END IF ! Approximative amplitude calculation

! ****************************************************************

      IF (.NOT.DYSEXPORT) THEN
       GOTO 100
      END IF

! Export part of the routine and exact calculation of amplitudes

! Read in the atomic overlap matrix, that will be needed below for
! for normalization of DOs
! (Code from mksxy)
      NSZZ=0
      DO 10 ISY=1,NSYM
        NO=NOSH(ISY)
        NB=NBASF(ISY)
        NSZZ=NSZZ+(NB*(NB+1))/2
        NSSQ=MAX(NSSQ,NB**2)
        NPROD=MAX(NPROD,NO*NB)
10    CONTINUE
      CALL GETMEM('SZZ   ','ALLO','REAL',LSZZ,NSZZ)
      IRC=-1
      IOPT=6
      ICMP=1
      ISYLAB=1
      CALL RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSZZ),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE SODYSORB ***'
        WRITE(6,*)'     OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF

! SZZ is originally given in triangular form, lets make it a full
! matrix for convenience
      NDUM=0
      DO ZJ=1,NZ
       DO ZI=1,ZJ
        SZZFULL(ZJ,ZI)=WORK(LSZZ+NDUM)
        SZZFULL(ZI,ZJ)=WORK(LSZZ+NDUM)
        NDUM=NDUM+1
       END DO
      END DO
      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NSZZ)

! ****************************************************************

C Multiply together with the SO eigenvector coefficients with the SF
C Dyson orbital coefficients in the atomic basis to obtain
C SO Dyson orbitals

      SODYSCIND=0 ! Orbital coeff. index
      ORBNUM=0 ! Dysorb index for given JSTATE
      SODYSCMOR=0.0D0
      SODYSCMOI=0.0D0
      DYSEN=0.0D0
      AMPS=0.0D0
      SODYSAMPS=0.0D0

      ! For all requested initial states J and all final states I
      DO JSTATE=1,DYSEXPSO
         DO ISTATE=JSTATE+1,NSS
      ! This loop nestling order is for some reason the reverse
      ! of what is done in GTDMCTL

        ! Reset values for next state combination
        DIJ=0.0D0
        DO NDUM=1,NZ
         SODYSCOFSR(NDUM)=0.0D0
         SODYSCOFSI(NDUM)=0.0D0
        END DO

        ! Iterate over the eigenvector components of both states
        DO JEIG=1,NSS

         ! Coefficient of first state
         INDJ=NSS*(JSTATE-1)+JEIG-1
         CJR=WORK(LUTOTR+INDJ)
         CJI=WORK(LUTOTI+INDJ)
         ! Find the corresponding SF states
         SFJ=IWORK(SO2SFNUM+JEIG-1)

         DO IEIG=1,NSS

          ! Coefficient of second state
          INDI=NSS*(ISTATE-1)+IEIG-1
          CIR=WORK(LUTOTR+INDI)
          CII=WORK(LUTOTI+INDI)
          ! Find the corresponding SF states
          SFI=IWORK(SO2SFNUM+IEIG-1)

          IF (DYSAMPS(SFJ,SFI).GT.1.0D-6) THEN
           ! Multiply together coefficients
           CREAL=CJR*CIR+CJI*CII
           CIMAG=CJR*CII-CJI*CIR
           ! Multiply with the corresponding SF Dyson orbital
           SODYSCOFSR=SODYSCOFSR+CREAL*SFDYS(:,SFJ,SFI)
           SODYSCOFSI=SODYSCOFSI+CIMAG*SFDYS(:,SFJ,SFI)
          END IF

         END DO ! IEIG
        END DO ! JEIG

! Normalize the overlap of SODYSCOFS expanded orbitals with the
! atomic overlap matrix SZZ to obtain correct amplitudes

          AMPLITUDE=0.0D0

          DO ZJ=1,NZ
           DO ZI=1,NZ
            AMPR=SODYSCOFSR(ZJ)*SODYSCOFSR(ZI)
     &            +SODYSCOFSI(ZJ)*SODYSCOFSI(ZI)
            AMPI=SODYSCOFSI(ZJ)*SODYSCOFSR(ZI)
     &            -SODYSCOFSR(ZJ)*SODYSCOFSI(ZI)
            AMPLITUDE=AMPLITUDE+(AMPR+AMPI)*SZZFULL(ZJ,ZI)
           END DO ! ZI
          END DO ! ZJ

          AMPLITUDE=SQRT(AMPLITUDE)
          SODYSAMPS(JSTATE,ISTATE)=AMPLITUDE
          SODYSAMPS(ISTATE,JSTATE)=AMPLITUDE

!+++  J. Creutzberg, J. Norell  - 2018 (.DysOrb and .molden export )
!     For each requested initial state JSTATE we will gather all the
!     obtained Dysorbs and export to a shared .DysOrb file and .molden
!     file.
!     Each file can however only contain NZ numbe of orbitals, so we
!     might have to split into several files IFILE
          IF ( ISTATE.EQ.(JSTATE+1) ) THEN
              IFILE=1
              SODYSCIND=0 ! Orbital coeff. index
              ORBNUM=0 ! Dysorb index for given JSTATE
              SODYSCMOR=0.0D0
              SODYSCMOI=0.0D0
              DYSEN=0.0D0
              AMPS=0.0D0
          END IF

          ! Export Re and Im part of the coefficients
          DO NDUM=1,NZ
             SODYSCIND=SODYSCIND+1
             SODYSCMOR(SODYSCIND)=SODYSCOFSR(NDUM)
             SODYSCMOI(SODYSCIND)=SODYSCOFSI(NDUM)
          END DO
          ORBNUM=ORBNUM+1
          DYSEN(ORBNUM)=SOENE(ISTATE)-SOENE(JSTATE)
          AMPS(ORBNUM)=AMPLITUDE*AMPLITUDE

! Write the Dysorbs from JSTATE to .DysOrb and .molden file
! (Enough to fill one file)
         IF(ORBNUM.EQ.NZ) THEN
          Call Dys_Interf(1,JSTATE,IFILE,NZ,SODYSCMOR,
     &        DYSEN,AMPS)
          Call Dys_Interf(2,JSTATE,IFILE,NZ,SODYSCMOI,
     &        DYSEN,AMPS)
          IFILE=IFILE+1
          SODYSCIND=0 ! Orbital coeff. index
          ORBNUM=0 ! Dysorb index for given JSTATE
          SODYSCMOR=0.0D0
          SODYSCMOI=0.0D0
          DYSEN=0.0D0
          AMPS=0.0D0
         END IF
! +++

        END DO ! ISTATE

! +++ J. Creutzberg, J. Norell - 2018
! Write the Dysorbs from JSTATE to .DysOrb and .molden file
! (All remaining, if any)
        IF(ORBNUM.GT.0) THEN
        Call Dys_Interf(1,JSTATE,IFILE,NZ,SODYSCMOR,
     &        DYSEN,AMPS)
        Call Dys_Interf(2,JSTATE,IFILE,NZ,SODYSCMOI,
     &        DYSEN,AMPS)
        END IF
! +++

       END DO ! JSTATE

100    CONTINUE

! ****************************************************************

C Free all the allocated memory

      CALL GETMEM('SO2SF','FREE','INTE',SO2SFNUM,NSS)

      RETURN
      END



