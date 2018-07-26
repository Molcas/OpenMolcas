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
      SUBROUTINE DO_SODYSORB(NSS,LUTOTR,LUTOTI,DYSAMPS,SFDYS,NZ,
     &                       SODYSAMPS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "Files.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='DO_SONATORB')
      LOGICAL   FIRSTSO
      INTEGER   SOTOT,SFTOT,SO2SFNUM,DUMMY
      INTEGER   NZ,LSZZ
      INTEGER   IDISK
      INTEGER   INDJ,INDI,SFI,SFJ,ZI,ZJ,NSZZ,NDUM
      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION SFDYS(NZ,NSTATE,NSTATE)
      DIMENSION SODYSAMPS(NSS,NSS)
      DIMENSION SODYSCOFS(NZ)
      DIMENSION SZZFULL(NZ,NZ)

C Calculates spin-orbit Dyson orbitals
C The routine was in large part adapted from DO_SONATORB

C 1. Perform a mapping from the SO eigenstates to the SF eigenstates
C 2. Pick up the SF Dyson orbitals in the atomic basis from disk
C 3. Build the SO Dyson orbitals from the SF ones
C 4. Save the SO Dyson orbitals to disk

      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*****************************************'
      WRITE(6,*) '* RUNNING SODYSORB CODE *****************'
      WRITE(6,*) '*****************************************'
      WRITE(6,*)

      SO2SFNUM=0
      CALL GETMEM('SO2SF','ALLO','INTE',SO2SFNUM,NSS)

C Setup SO2SFNUM list which contains the original SF state numbers
C as a function of the SO state number
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

C Write out SO2SFNUM for debugging
!      WRITE(*,*)'---------------------------'
!      WRITE(*,*)'------- SO2SFNUM ----------'
!      WRITE(*,*)'---------------------------'
!      DO ISTATE=1,NSS
!        WRITE(*,'(F5.2)',advance="no"),WORK(SO2SFNUM+ISTATE-1)
!       WRITE(*,*)
!      END DO

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

C Write out V-matrix for debugging
!      WRITE(*,*)'---------------------------'
!      WRITE(*,*)'------- V-matrix ----------'
!      WRITE(*,*)'---------------------------'
!      DO JSTATE=1,NSS
!       DO ISTATE=1,NSS
!        WRITE(*,'(F5.2)',advance="no"),WORK(LVMAT-1+NSS*(ISTATE-1)
!     &    +JSTATE)
!       END DO
!       WRITE(*,*)
!      END DO

! ****************************************************************

C Now multiply V with the SO eigenvector matrices to obtain the
C the SO eigenstates in the (multiplicity expanded) spin-free basis.
      CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,WORK(LUTOTR),NSS,0.0d0,
     &      WORK(LUMATR),NSS)
      CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,WORK(LUTOTI),NSS,0.0d0,
     &      WORK(LUMATI),NSS)

C Write out V*E-matrix for debugging
!      WRITE(*,*)'---------------------------'
!      WRITE(*,*)'------- V*E-matrix --------'
!      WRITE(*,*)'---------------------------'
!      DO ISTATE=1,NSS
!       DO JSTATE=1,NSS
!        NINDEX=NSS*(ISTATE-1)+JSTATE-1
!        WRITE(*,'(F5.2,A,F4.2,A)',advance="no")
!     &        WORK(LUTOTR+NINDEX),'+',WORK(LUOTI+NINDEX),'i'
!       END DO
!       WRITE(*,*)
!      END DO

! ****************************************************************

C Now read in all the previously saved SF Dyson orbitals in the
C atomic basis from disk

      DO JSTATE=1,NSTATE
       DO ISTATE=JSTATE,NSTATE
       IF (DYSAMPS(ISTATE,JSTATE).GT.1.0D-6) THEN
        IDISK=IWORK(LIDDYS+(ISTATE-1)*NSTATE+JSTATE-1)
        CALL DDAFILE(LUDYS,2,SFDYS(:,ISTATE,JSTATE),NZ,IDISK)
       END IF
       END DO
      END DO

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

C Write out SZZ for debugging
!      WRITE(*,*)'SZZ:'
!      DO ZJ=1,NZ
!       DO ZI=1,NZ
!        WRITE(*,'(F6.2)')SZZFULL(ZJ,ZI)
!       END DO
!       WRITE(*,*)
!      END DO

! ****************************************************************

C Multiply together with the SO eigenvector coefficients to get the
C SO coefficients

      ! For all possible SO state combinations
      DO JSTATE=1,NSS
       DO ISTATE=1,JSTATE

        ! Reset values for next state combination
        DIJ=0.0D0
        DO NDUM=1,NZ
         SODYSCOFS(NDUM)=0.0D0
        END DO

        ! Iterate over the eigenvector components of both states
        DO JEIG=1,NSS
         DO IEIG=1,NSS

          ! Lets start with real values for testing and add imag part
          ! later
          ! Coefficient of first state
          INDJ=NSS*(JSTATE-1)+JEIG-1
          CJ=WORK(LUTOTR+INDJ)
          ! Coefficient of second state
          INDI=NSS*(ISTATE-1)+IEIG-1
          CI=WORK(LUTOTR+INDI)
          ! Find the corresponding SF states
          SFJ=WORK(SO2SFNUM+JEIG-1)
          SFI=WORK(SO2SFNUM+IEIG-1)
          ! Multiply coefficients with the corresponding SF Dyson orbital
          SODYSCOFS=SODYSCOFS+(CJ*CI)*SFDYS(:,SFJ,SFI)

         END DO
        END DO

C Write out SODYSCOFS for debugging
!        WRITE(*,*)
!        WRITE(*,*)'SODYSCOFS='
!        DO NDUM=1,NZ
!         WRITE(*,'(F6.2)',advance='no')SODYSCOFS(NDUM)
!        END DO

! Normalize the overlap of SODYSCOFS expanded orbitals with the
! atomic overlap matrix SZZ to obtain correct amplitudes

          AMPLITUDE=0.0D0

          DO ZJ=1,NZ
           DO ZI=1,NZ
            ! SZZ normalization must be added!
            AMPLITUDE=AMPLITUDE+SODYSCOFS(ZJ)*SODYSCOFS(ZI)
     &                *SZZFULL(ZJ,ZI)
           END DO
          END DO

          AMPLITUDE=SQRT(AMPLITUDE)
          SODYSAMPS(JSTATE,ISTATE)=AMPLITUDE
          SODYSAMPS(ISTATE,JSTATE)=AMPLITUDE

        END DO
       END DO

! The coefficients could be used for some kind of analysis,
! but for now we will not save them before going to the next
! state combination

C Write out SODYSAMPS for debugging
!      WRITE(*,*)'SODYSAMPS='
!      DO ISS=1,NSS
!       DO JSS=1,NSS
!        WRITE(*,*)SODYSAMPS(JSS,ISS)
!       END DO
!       WRITE(*,*)
!      END DO

! ****************************************************************

C Free all the allocated memory

      CALL GETMEM('SO2SF','FREE','INTE',SO2SFNUM,NSS)

      CALL GETMEM('UMATR2','FREE','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','FREE','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','FREE','REAL',LVMAT,NSS**2)


      RETURN
      END



