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
      SUBROUTINE DO_SONATORB(NSS, USOR, USOI)
      use rassi_aux, only: ipglob
      use rassi_global_arrays, only: JBNUM, EIGVEC
      use cntrl, only: SONAT, SONATNSTATE,
     &                      SODIAG, SODIAGNSTATE
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, NOSO, IfCurd, MLTPLT
      use rassi_data, only: NBTRI
      IMPLICIT None
      INTEGER NSS
      Real*8 USOR(NSS,NSS), USOI(NSS,NSS)
      Real*8 IDENTMAT(3,3)
      Real*8, Allocatable:: UMATR(:), UMATI(:), VMAT(:,:)
      Real*8, Allocatable:: DMATTMP(:)
      Integer ISS, ISTATE, JOB1, MPLET1, MSPROJ1,
     &        JSS, JSTATE, JOB2, MPLET2, MSPROJ2,
     &        I, INATSTATE, IOPT, IC
      REAL*8 DUM1, DUM2, DUM3, DUM4, DUM5, DUM6

c Calculates natural orbitals, including spinorbit effects
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*****************************************'
      WRITE(6,*) '* RUNNING SONATORB CODE *****************'
      WRITE(6,*) '*****************************************'
      WRITE(6,*)

      IDENTMAT(:,:)=0.0D0
      IDENTMAT(1,1)=1.0D0
      IDENTMAT(2,2)=1.0D0
      IDENTMAT(3,3)=1.0D0

      CALL mma_allocate(UMATR,NSS**2,Label='UMATR')
      CALL mma_allocate(UMATI,NSS**2,Label='UMATI')
      CALL mma_allocate(VMAT,NSS,NSS,Label='VMAT')

      VMAT(:,:)=0.0D0

c transform V matrix in SF basis to spin basis
c This was taken from smmat.f and modified slightly
      ISS=0
      DO ISTATE=1,NSTATE
       JOB1=JBNUM(ISTATE)
       MPLET1=MLTPLT(JOB1)

       DO MSPROJ1=-MPLET1+1,MPLET1-1,2
        ISS=ISS+1
        JSS=0

        DO JSTATE=1,NSTATE
         JOB2=JBNUM(JSTATE)
         MPLET2=MLTPLT(JOB2)

         DO MSPROJ2=-MPLET2+1,MPLET2-1,2
          JSS=JSS+1

          IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2) THEN
           VMAT(ISS,JSS)=EIGVEC(JSTATE,ISTATE)
          END IF ! IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2)
         END DO ! DO MSPROJ2=-MPLET2+1,MPLET2-1,2
        END DO ! end DO JSTATE=1,NSTATE
       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE



c combine this matrix with the SO eigenvector matrices
      IF(.not.NOSO) THEN
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,VMAT,NSS,USOR,NSS,0.0d0,
     &      UMATR,NSS)
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,VMAT,NSS,USOI,NSS,0.0d0,
     &      UMATI,NSS)
      ELSE
c Spinorbit contributions to this are disabled
        CALL DCOPY_(NSS,VMAT,1,UMATR,1)
        CALL DCOPY_(NSS,[0.0d0],0,UMATI,1)
      END IF

c Holds the density matrices for all three directions
      CALL mma_allocate(DMATTMP,6*NBTRI,Label='DMATTMP')

c SONATNSTATE = number of states to calculate.
c These states are stored beginning in SONAT
      DO I=1,SONATNSTATE
        INATSTATE=SONAT(I)

        WRITE(6,*)
        WRITE(6,*) "CALCULATING NAT ORBITALS FOR SSTATE: ",INATSTATE
        IF(INATSTATE.GT.NSS.OR.INATSTATE.LE.0) THEN
          WRITE(6,*) "...WHICH DOES NOT EXIST!"
          CALL ABEND()
        END IF
        WRITE(6,*)

C Calculate overall density, store in DMATTMP
        iOpt=0
        CALL SONATORBM('HERMSING',UMATR,UMATI,
     &                 INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                 DMATTMP)


C Integrate for the expectation value
        IF(IPGLOB.ge.3) THEN
          IC=1
          iOpt=0
          CALL SONATORBM_INT(DMATTMP,'MLTPL  0',IC,'HERMSING',
     &                      INATSTATE,INATSTATE,iOpt,IDENTMAT,
     &                      DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)

C          CALL ADD_INFO('MLTPL0SING_INT3',DUM3,1,6)
C          CALL ADD_INFO('MLTPL0SING_INT6',DUM6,1,6)
        END IF

C Create SONATTDENS total density orbital file for this (I,I) state
        CALL SONATORB_PLOT(DMATTMP,'SONATTDENS','HERMSING',
     &                     INATSTATE,INATSTATE)

        IF(IPGLOB.ge.4) THEN
          CALL SONATORB_CPLOT(DMATTMP,'TDENSTESTX','HERMSING',
     &                       INATSTATE,INATSTATE)
        END IF




C Calculate spin density, store in LDMATTMP
        iOpt=0
        CALL SONATORBM('HERMTRIP',UMATR,UMATI,
     &                 INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                 DMATTMP)


C Integrate for the expectation value
        IF(IPGLOB.ge.3) THEN
          IC=1
          iOpt=0
          CALL SONATORBM_INT(DMATTMP,'MLTPL  0',IC,'HERMTRIP',
     &                      INATSTATE,INATSTATE,iOpt,IDENTMAT,
     &                      DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
C          CALL ADD_INFO('MLTPL0TRIP_INT3',DUM3,1,6)
C          CALL ADD_INFO('MLTPL0TRIP_INT6',DUM6,1,6)
        END IF


C Create SONATSDENS spin density orbital file for this (I,I) state
        CALL SONATORB_PLOT(DMATTMP,'SONATSDENS','HERMTRIP',
     &                     INATSTATE,INATSTATE)

        IF(IPGLOB.ge.4) THEN
          CALL SONATORB_CPLOT(DMATTMP,'SDENSTESTX','HERMTRIP',
     &                       INATSTATE,INATSTATE)
        END IF


c Type 2 - current density
        IF(IFCURD) THEN
          iOpt=0
          CALL SONATORBM('ANTISING',UMATR,UMATI,
     &                   INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                   DMATTMP)


          IF(IPGLOB.ge.3) THEN
            IC=1
            iOpt=0
            CALL SONATORBM_INT(DMATTMP,'ANGMOM  ',IC,'ANTISING',
     &                        INATSTATE,INATSTATE,iOpt,IDENTMAT,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
C            CALL ADD_INFO('CURD1_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD1_INT6',DUM6,1,6)

            IC=2
            iOpt=0
            CALL SONATORBM_INT(DMATTMP,'ANGMOM  ',IC,'ANTISING',
     &                        INATSTATE,INATSTATE,iOpt,IDENTMAT,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)

C            CALL ADD_INFO('CURD2_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD2_INT6',DUM6,1,6)

            IC=3
            iOpt=0
            CALL SONATORBM_INT(DMATTMP,'ANGMOM  ',IC,'ANTISING',
     &                        INATSTATE,INATSTATE,iOpt,IDENTMAT,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
C            CALL ADD_INFO('CURD3_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD3_INT6',DUM6,1,6)

          END IF


          CALL SONATORB_CPLOT(DMATTMP,'SONATLDENS','ANTISING',
     &                        INATSTATE,INATSTATE)
        END IF


      END DO

      CALL mma_deallocate(DMATTMP)
      CALL mma_deallocate(SONAT)


c perform the state diagonalization similar to
c what is done in single_aniso
      IF(SODIAGNSTATE.GT.0) THEN

c This actually does all the work
        CALL mkSODIAG(UMATR, UMATI, NSS)

c This is only allocated if SODIAGNSTATE.GT.0
        CALL mma_deallocate(SODIAG)
      END IF

      CALL mma_deallocate(UMATR)
      CALL mma_deallocate(UMATI)
      CALL mma_deallocate(VMAT)

      END SUBROUTINE DO_SONATORB
