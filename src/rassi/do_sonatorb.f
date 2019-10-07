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
      SUBROUTINE DO_SONATORB(NSS, LUTOTR, LUTOTI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='DO_SONATORB')
      Real*8 IDENTMAT(3,3)

c Calculates natural orbitals, including spinorbit effects
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*****************************************'
      WRITE(6,*) '* RUNNING SONATORB CODE *****************'
      WRITE(6,*) '*****************************************'
      WRITE(6,*)

      IDENTMAT(:,:)=0.0D0
      FOR ALL (I=1:3) IDENTMAT(I,I)=1.0D0

      CALL GETMEM('UMATR2','ALLO','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','ALLO','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','ALLO','REAL',LVMAT,NSS**2)

      CALL DCOPY_(NSS**2,[0.0d0],0,WORK(LVMAT),1)

c transform V matrix in SF basis to spin basis
c This was taken from smmat.f and modified slightly
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



c combine this matrix with the SO eigenvector matrices
      IF(.not.NOSO) THEN
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,WORK(LUTOTR),NSS,0.0d0,
     &      WORK(LUMATR),NSS)
        CALL DGEMM_('N','N',NSS,NSS,NSS,
     &      1.0d0,WORK(LVMAT),NSS,WORK(LUTOTI),NSS,0.0d0,
     &      WORK(LUMATI),NSS)
      ELSE
c Spinorbit contributions to this are disabled
        CALL DCOPY_(NSS,WORK(LVMAT),1,WORK(LUMATR),1)
        CALL DCOPY_(NSS,[0.0d0],0,WORK(LUMATI),1)
      END IF

c Holds the density matrices for all three directions
      CALL GETMEM('DMATTMP','ALLO','REAL',LDMATTMP,6*NBTRI)

c SONATNSTATE = number of states to calculate.
c These states are stored beginning in IWORK(LSONAT)
      DO I=1,SONATNSTATE
        INATSTATE=IWORK(LSONAT-1+I)

        WRITE(6,*)
        WRITE(6,*) "CALCULATING NAT ORBITALS FOR SSTATE: ",INATSTATE
        IF(INATSTATE.GT.NSS.OR.INATSTATE.LE.0) THEN
          WRITE(6,*) "...WHICH DOES NOT EXIST!"
          CALL ABEND()
        END IF
        WRITE(6,*)

C Calculate overall density, store in WORK(LDMATTMP)
        CALL SONATORB('HERMSING',WORK(LUMATR),WORK(LUMATI),
     &                INATSTATE,INATSTATE,NSS,WORK(LDMATTMP))


C Integrate for the expectation value
        IF(IPGLOB.ge.VERBOSE) THEN
          IC=1
          iOpt=0
          CALL SONATORBM_INT(WORK(LDMATTMP),'MLTPL  0',IC,'HERMSING',
     &                      INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                      DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)

C          CALL ADD_INFO('MLTPL0SING_INT3',DUM3,1,6)
C          CALL ADD_INFO('MLTPL0SING_INT6',DUM6,1,6)
        END IF

C Create SONATTDENS total density orbital file for this (I,I) state
        CALL SONATORB_PLOT(WORK(LDMATTMP),'SONATTDENS','HERMSING',
     &                     INATSTATE,INATSTATE)

        IF(IPGLOB.ge.DEBUG) THEN
          CALL SONATORB_CPLOT(WORK(LDMATTMP),'TDENSTESTX','HERMSING',
     &                       INATSTATE,INATSTATE)
        END IF




C Calculate spin density, store in LDMATTMP
        CALL SONATORB('HERMTRIP',WORK(LUMATR),WORK(LUMATI),
     &                INATSTATE,INATSTATE,NSS,WORK(LDMATTMP))


C Integrate for the expectation value
        IF(IPGLOB.ge.VERBOSE) THEN
          IC=1
          iOpt=0
          CALL SONATORBM_INT(WORK(LDMATTMP),'MLTPL  0',IC,'HERMTRIP',
     &                      INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                      DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
C          CALL ADD_INFO('MLTPL0TRIP_INT3',DUM3,1,6)
C          CALL ADD_INFO('MLTPL0TRIP_INT6',DUM6,1,6)
        END IF


C Create SONATSDENS spin density orbital file for this (I,I) state
        CALL SONATORB_PLOT(WORK(LDMATTMP),'SONATSDENS','HERMTRIP',
     &                     INATSTATE,INATSTATE)

        IF(IPGLOB.ge.DEBUG) THEN
          CALL SONATORB_CPLOT(WORK(LDMATTMP),'SDENSTESTX','HERMTRIP',
     &                       INATSTATE,INATSTATE)
        END IF


c Type 2 - current density
        IF(IFCURD) THEN
          CALL SONATORB('ANTISING',WORK(LUMATR),WORK(LUMATI),
     &                   INATSTATE,INATSTATE,NSS,WORK(LDMATTMP))


          IF(IPGLOB.ge.VERBOSE) THEN
            IC=1
            iOpt=0
            CALL SONATORBM_INT(WORK(LDMATTMP),'ANGMOM  ',IC,'ANTISING',
     &                        INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
C            CALL ADD_INFO('CURD1_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD1_INT6',DUM6,1,6)

            IC=2
            iOpt=0
            CALL SONATORBM_INT(WORK(LDMATTMP),'ANGMOM  ',IC,'ANTISING',
     &                        INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)

C            CALL ADD_INFO('CURD2_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD2_INT6',DUM6,1,6)

            IC=3
            iOpt=0
            CALL SONATORBM_INT(WORK(LDMATTMP),'ANGMOM  ',IC,'ANTISING',
     &                        INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
C            CALL ADD_INFO('CURD3_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD3_INT6',DUM6,1,6)

          END IF


          CALL SONATORB_CPLOT(WORK(LDMATTMP),'SONATLDENS','ANTISING',
     &                        INATSTATE,INATSTATE)
        END IF


      END DO

      CALL GETMEM('DMATTMP','FREE','REAL',LDMATTMP,6*NBTRI)
      CALL GETMEM('SONATS','FREE','INTE',LSONAT,SONATNSTATE)


c perform the state diagonalization similar to
c what is done in single_aniso
      IF(SODIAGNSTATE.GT.0) THEN

!       Triangular part of a matrix
        LSOSIZ=SODIAGNSTATE*(SODIAGNSTATE+1)

c This actually does all the work
        CALL SODIAG(WORK(LUMATR), WORK(LUMATI), NSS)

c This is only allocated if SODIAGNSTATE.GT.0
        CALL GETMEM('SODIAG','FREE','INTE',LSODIAG,SODIAGNSTATE)
      END IF

      CALL GETMEM('UMATR2','FREE','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','FREE','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','FREE','REAL',LVMAT,NSS**2)

      RETURN
      END
