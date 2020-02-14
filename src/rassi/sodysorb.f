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
* Copyright (C) 1987, Per Ake Malmqvist                                *
*               2018, Jesper Norell                                    *
*               2018, Joel Creutzberg                                  *
************************************************************************
      SUBROUTINE SODYSORB(NSS,USOR,USOI,DYSAMPS,NZ,SOENE)
      use rassi_global_arrays, only: SFDYS, SODYSAMPS,
     &                               SODYSAMPSR, SODYSAMPSI, JBNUM

      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "Files.fh"

      REAL*8 USOR(NSS,NSS), USOI(NSS,NSS)
      ! Arrays, bounds, and indices
      INTEGER   MSPROJS
      REAL*8    MSPROJI,MSPROJJ
      INTEGER   SOTOT,SFTOT,SO2SFNUM
      INTEGER   NZ,LSZZ,ORBNUM
      INTEGER   SODYSCIND
      INTEGER   SFI,SFJ,ZI,ZJ,NSZZ,NDUM

      ! Arrays for calculation of amplitudes
      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION SODYSCOFSR(NZ),SODYSCOFSI(NZ)
      DIMENSION SZZFULL(NZ,NZ)

      ! Arrays for orbital export
      DIMENSION SODYSCMOR(NZ*NSS)
      DIMENSION SODYSCMOI(NZ*NSS)
      DIMENSION SOENE(NSS)
      DIMENSION DYSEN(NSS)
      DIMENSION AMPS(NSS)
      Character*30 Filename

! +++ J.Norell 2018

C Calculates spin-orbit Dyson orbitals
C The routine was in some part adapted from DO_SONATORB

C Computes SO Dyson amplitudes by expanding the SF results with
C the SO eigenvectors (of the complex Hamiltonian)
C 1. (Fast): Compute the SO amplitudes directly from the SF amplitudes
C    (approximation) for all states
C 2. (Slower): Compute the full SO Dyson orbitals for the requested
C    initial states and export them to .molden format. SO amplitudes
C    are correctly calculated for these states.

! ****************************************************************

C Setup SO2SFNUM list which contains the original SF state numbers
C as a function of the SO state number
C And MSPROJS which saves their ms projections for later use

      SO2SFNUM=0
      CALL GETMEM('SO2SF','ALLO','INTE',SO2SFNUM,NSS)
      MSPROJS=0
      CALL GETMEM('MSPROJS','ALLO','REAL',MSPROJS,NSS)
      ISS=0
      SOTOT=0
      SFTOT=0
      DO ISTATE=1,NSTATE
       JOB1=JBNUM(ISTATE)
       MPLET1=MLTPLT(JOB1)
       SFTOT=SFTOT+1

       DO MSPROJ=-MPLET1+1,MPLET1-1,2
        SOTOT=SOTOT+1
        IWORK(SO2SFNUM+SOTOT-1)=SFTOT
        WORK(MSPROJS+SOTOT-1)=MSPROJ

       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE

      IF(.NOT.DYSEXPORT) THEN ! Approximative amplitude calculation

C Now construct the SF dysamps in the multiplicity expanded basis
C (initially all real, therefore put into SODYSAMPSR)
      SODYSAMPSR(:,:)=0.0D0
      SODYSAMPSI(:,:)=0.0D0
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
      CALL ZTRNSF(NSS,USOR,USOI,SODYSAMPSR,SODYSAMPSI)

C Compute the magnitude of the complex amplitudes as an approximation
      SODYSAMPSR(:,:)=SODYSAMPSR*SODYSAMPSR
      SODYSAMPSI(:,:)=SODYSAMPSI*SODYSAMPSI
      SODYSAMPS(:,:)=SQRT(SODYSAMPSR+SODYSAMPSI)

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
      NSSQ=0
      NPROD=0
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

! SZZ is originally given in symmetry-blocked triangular form,
! lets make it a full matrix for convenience
      SZZFULL=0.0D0
      NDUM=0
      NOFF=0
      DO ISY=1,NSYM
       NB=NBASF(ISY)
       DO ZJ=1,NB
        DO ZI=1,ZJ
         SZZFULL(ZJ+NOFF,ZI+NOFF)=WORK(LSZZ+NDUM)
         SZZFULL(ZI+NOFF,ZJ+NOFF)=WORK(LSZZ+NDUM)
         NDUM=NDUM+1
        END DO
       END DO
       NOFF=NOFF+NB
      END DO
      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NSZZ)

! ****************************************************************

C Multiply together with the SO eigenvector coefficients with the SF
C Dyson orbital coefficients in the atomic basis to obtain the full
C SO Dyson orbitals

C Multiply together with the SO eigenvector coefficients with the SF
C Dyson orbital coefficients in the atomic basis to obtain
C SO Dyson orbitals

      SODYSAMPS(:,:)=0.0D0
      ! For all requested initial states J and all final states I
      DO JSTATE=1,DYSEXPSO

         ! For each initial state JSTATE up to DYSEXPSFSO we will
         ! gather all the obtained Dysorbs
         ! and export to a shared .molden file
         IFILE=1
         SODYSCIND=0 ! Orbital coeff. index
         ORBNUM=0 ! Dysorb index for given JSTATE
         SODYSCMOR=0.0D0 ! Real orbital coefficients
         SODYSCMOI=0.0D0 ! Imaginary orbital coefficients
         DYSEN=0.0D0 ! Orbital energies
         AMPS=0.0D0 ! Transition amplitudes (shown as occupations)

         DO ISTATE=JSTATE+1,NSS

          ! Reset values for next state combination
          SODYSCOFSR=0.0D0
          SODYSCOFSI=0.0D0

          ! Iterate over the eigenvector components of both states
          DO JEIG=1,NSS

           ! Coefficient of first state
           CJR=USOR(JEIG,JSTATE)
           CJI=USOI(JEIG,JSTATE)
           ! Find the corresponding SF states
           SFJ=IWORK(SO2SFNUM+JEIG-1)

           DO IEIG=1,NSS

            ! Coefficient of second state
            CIR=USOR(IEIG,ISTATE)
            CII=USOI(IEIG,ISTATE)
            ! Find the corresponding SF states
            SFI=IWORK(SO2SFNUM+IEIG-1)

            ! Check change in ms projection
            MSPROJJ=WORK(MSPROJS+JEIG-1)
            MSPROJI=WORK(MSPROJS+IEIG-1)
            ! Check |delta ms|=0.5 selection rule
            IF(ABS(MSPROJJ-MSPROJI).NE.1) THEN
             CYCLE
            END IF

            IF (DYSAMPS(SFJ,SFI).GT.1.0D-5) THEN
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

          ! Export Re and Im part of the coefficients
          IF (AMPLITUDE.GT.1.0D-5) THEN
           DO NDUM=1,NZ
              SODYSCIND=SODYSCIND+1
              SODYSCMOR(SODYSCIND)=SODYSCOFSR(NDUM)
              SODYSCMOI(SODYSCIND)=SODYSCOFSI(NDUM)
           END DO
           ORBNUM=ORBNUM+1
           DYSEN(ORBNUM)=SOENE(ISTATE)-SOENE(JSTATE)
           AMPS(ORBNUM)=AMPLITUDE*AMPLITUDE
          END IF

        END DO ! ISTATE

! If at least one orbital was found, export it/them
        IF(ORBNUM.GT.0) THEN
         Write(filename,'(A16,I0,A3)') 'Dyson.SO.molden.',JSTATE,'.Re'
         Call Molden_DysOrb(0,filename,DYSEN,AMPS,SODYSCMOR,ORBNUM,NZ)
         Write(filename,'(A16,I0,A3)') 'Dyson.SO.molden.',JSTATE,'.Im'
         Call Molden_DysOrb(0,filename,DYSEN,AMPS,SODYSCMOI,ORBNUM,NZ)
        END IF

       END DO ! JSTATE

100    CONTINUE

! ****************************************************************

C Free all the allocated memory

      CALL GETMEM('SO2SF','FREE','INTE',SO2SFNUM,NSS)
      CALL GETMEM('MSPROJS','FREE','REAL',MSPROJS,NSS)

      RETURN
      END
