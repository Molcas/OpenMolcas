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
     &                       SODYSAMPS,SOENE)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "Files.fh"
      CHARACTER*16 ROUTINE
      LOGICAL   FIRSTSO
      INTEGER   SOTOT,SFTOT,SO2SFNUM,DUMMY
      INTEGER   NZ,LSZZ,ORBNUM
      INTEGER   IDISK
      INTEGER   SODYSCIND,ENIND
      INTEGER   INDJ,INDI,SFI,SFJ,ZI,ZJ,NSZZ,NDUM

      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION SFDYS(NZ,NSTATE,NSTATE)
      DIMENSION SODYSAMPS(NSS,NSS)
      DIMENSION SODYSCOFSR(NZ),SODYSCOFSI(NZ)
      DIMENSION SZZFULL(NZ,NZ)

      ! Joel
      DIMENSION SODYSCMO(NZ*NZ)
      DIMENSION SOENE(NZ)
      DIMENSION DYSEN(NZ)
      DIMENSION AMPS(NZ)

! +++ J.Norell 2018

C Calculates spin-orbit Dyson orbitals
C The routine was in some part adapted from DO_SONATORB

C 1. Perform a mapping from the SO eigenstates to the SF eigenstates
C 2. Pick up the SF Dyson orbitals in the atomic basis from disk
C 3. Build the SO Dyson orbitals (and their amplitudes) from the SF ones
C 4. (Optional) export the SO orbitals

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
        WORK(SO2SFNUM+SOTOT-1)=SFTOT

       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE

! **** This part seems unecessary, as it essentially appears *******
! **** to create an identity matrix, but lets leave it for now!

      CALL GETMEM('UMATR2','ALLO','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','ALLO','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','ALLO','REAL',LVMAT,NSS**2)
      CALL DCOPY_(NSS**2,0.0d0,0,WORK(LVMAT),1)

C Setup transformation matrix V that expresses the initial SO
C states in the (multiplicity expanded) SF basis
      ISS=0
      DO ISTATE=1,NSTATE
       JOB1=iWork(lJBNUM+ISTATE-1)
       MPLET1=MLTPLT(JOB1)
       S1=0.5D0*DBLE(MPLET1-1)

       DO MSPROJ1=-MPLET1+1,MPLET1-1,2
        SM1=0.5D0*DBLE(MSPROJ1)
        ISS=ISS+1
        JSS=0

        DO JSTATE=1,NSTATE
         JOB2=iWork(lJBNUM+JSTATE-1)
         MPLET2=MLTPLT(JOB2)
         S2=0.5D0*DBLE(MPLET2-1)

         DO MSPROJ2=-MPLET2+1,MPLET2-1,2
          SM2=0.5D0*DBLE(MSPROJ2)
          JSS=JSS+1

          IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2) THEN
           IJ=(JSS-1)*NSS+ISS
           WORK(LVMAT-1+IJ)=WORK(LEIGVEC+(ISTATE-1)*NSTATE+JSTATE-1)
          END IF ! IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2)
         END DO ! DO MSPROJ2=-MPLET2+1,MPLET2-1,2
        END DO ! end DO JSTATE=1,NSTATE
       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE

! ****************************************************************

C Now multiply V with the SO eigenvector matrices to obtain the
C the SO eigenstates in the (multiplicity expanded) spin-free basis.
      CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,WORK(LUTOTR),NSS,0.0d0,
     &      WORK(LUMATR),NSS)
      CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,WORK(LUTOTI),NSS,0.0d0,
     &      WORK(LUMATI),NSS)

! ****************************************************************

C Now read in all the previously saved SF Dyson orbitals in the
C atomic basis from disk

      DO JSTATE=1,NSTATE
       DO ISTATE=JSTATE+1,NSTATE
       ! This loop nestling order is for some reason the reverse
       ! of what is done in GTDMCTL
        IF (DYSAMPS(JSTATE,ISTATE).GT.1.0D-6) THEN
         IDISK=IWORK(LIDDYS+(ISTATE-1)*NSTATE+JSTATE-1)
         CALL DDAFILE(LUDYS,2,SFDYS(:,JSTATE,ISTATE),NZ,IDISK)
         ! Loops over states are performed triangularly, but
         ! permutation of degenerate states in the SO part
         ! might 'escape' this, therefore we fill out the
         ! full matrix to be safe.
!         SFDYS(:,ISTATE,JSTATE)=SFDYS(:,JSTATE,ISTATE)
        ELSE
         DO NDUM=1,NZ
          SFDYS(NDUM,JSTATE,ISTATE)=0.0D0
          SFDYS(NDUM,ISTATE,JSTATE)=0.0D0
         END DO ! NDUM
        END IF ! AMP Threshold
       END DO ! ISTATE
      END DO ! JSTATE

! ****************************************************************

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

! ****************************************************************

C Multiply together with the SO eigenvector coefficients to get the
C SO coefficients

      SODYSCIND=0 ! Orbital coeff. index
      ORBNUM=0 ! Dysorb index for given JSTATE
      SODYSCMO=0.0D0
      DYSEN=0.0D0
      AMPS=0.0D0

      ! For all possible SO state combinations
      DO JSTATE=1,NSS
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
         DO IEIG=1,NSS

          ! Coefficient of first state
          INDJ=NSS*(JSTATE-1)+JEIG-1
          CJR=WORK(LUTOTR+INDJ)
          CJI=WORK(LUTOTI+INDJ)
          ! Coefficient of second state
          INDI=NSS*(ISTATE-1)+IEIG-1
          CIR=WORK(LUTOTR+INDI)
          CII=WORK(LUTOTI+INDI)
          ! Find the corresponding SF states
          SFJ=WORK(SO2SFNUM+JEIG-1)
          SFI=WORK(SO2SFNUM+IEIG-1)
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
!     For each initial state JSTATE we will gather all the obtained Dysorbs
!     and export to a shared .DysOrb file and .molden file if
!     requested
          IF (.NOT.DYSEXPORT) THEN
              CYCLE
          ELSE IF (JSTATE.GT.DYSEXPSO) THEN
              CYCLE
          END IF
          IF ( ISTATE.EQ.(JSTATE+1) ) THEN
              SODYSCIND=0 ! Orbital coeff. index
              ORBNUM=0 ! Dysorb index for given JSTATE
              SODYSCMO=0.0D0
              DYSEN=0.0D0
              AMPS=0.0D0
          END IF

         IF (AMPLITUDE.GT.1.0D-6) THEN
          ! For now we are only exporting the real part
          ! of the orbitals, this should be updated later
          DO NDUM=1,NZ
             SODYSCIND=SODYSCIND+1
             SODYSCMO(SODYSCIND)=SODYSCOFSR(NDUM)
          END DO
          ORBNUM=ORBNUM+1
          DYSEN(ORBNUM)=SOENE(ISTATE)-SOENE(JSTATE)
          AMPS(ORBNUM)=AMPLITUDE*AMPLITUDE
         END IF
! +++

        END DO ! ISTATE

! +++ J. Creutzberg, J. Norell - 2018
! Write the Dysorbs from JSTATE to .DysOrb and .molden file
         IF (.NOT.DYSEXPORT) THEN
             CYCLE
         ELSE IF (JSTATE.GT.DYSEXPSO) THEN
             CYCLE
         END IF
         Call Dys_Interf(.TRUE.,JSTATE,NZ,SODYSCMO,
     &        DYSEN,AMPS)
! +++

       END DO ! JSTATE

! ****************************************************************

C Free all the allocated memory

      CALL GETMEM('SO2SF','FREE','INTE',SO2SFNUM,NSS)

      CALL GETMEM('UMATR2','FREE','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','FREE','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','FREE','REAL',LVMAT,NSS**2)
      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NSZZ)


      RETURN
      END



