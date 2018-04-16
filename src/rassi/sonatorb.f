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

c Calculates natural orbitals, including spinorbit effects
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*****************************************'
      WRITE(6,*) '* RUNNING SONATORB CODE *****************'
      WRITE(6,*) '*****************************************'
      WRITE(6,*)



      CALL GETMEM('UMATR2','ALLO','REAL',LUMATR,NSS**2)
      CALL GETMEM('UMATI2','ALLO','REAL',LUMATI,NSS**2)
      CALL GETMEM('EIGVEC2','ALLO','REAL',LVMAT,NSS**2)

      CALL DCOPY_(NSS**2,0.0d0,0,WORK(LVMAT),1)

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
        CALL DCOPY_(NSS,0.0d0,0,WORK(LUMATI),1)
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
          CALL SONATORB_INT(WORK(LDMATTMP),'MLTPL  0',1,'HERMSING',
     &                      INATSTATE,INATSTATE,NSS,
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
          CALL SONATORB_INT(WORK(LDMATTMP),'MLTPL  0',1,'HERMTRIP',
     &                      INATSTATE,INATSTATE,NSS,
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
            CALL SONATORB_INT(WORK(LDMATTMP),'ANGMOM  ',1,'ANTISING',
     &                        INATSTATE,INATSTATE,NSS,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
C            CALL ADD_INFO('CURD1_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD1_INT6',DUM6,1,6)

            CALL SONATORB_INT(WORK(LDMATTMP),'ANGMOM  ',2,'ANTISING',
     &                        INATSTATE,INATSTATE,NSS,
     &                        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)

C            CALL ADD_INFO('CURD2_INT3',DUM3,1,6)
C            CALL ADD_INFO('CURD2_INT6',DUM6,1,6)

            CALL SONATORB_INT(WORK(LDMATTMP),'ANGMOM  ',3,'ANTISING',
     &                        INATSTATE,INATSTATE,NSS,
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

        LSOSIZ=SODIAGNSTATE*(SODIAGNSTATE+1) ! Triangular part of a matrix

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



      SUBROUTINE SONATORB(CHARTYPE,
     &                   USOR,USOI,ASS,BSS,NSS,
     &                   DENSOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SONATORB')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION USOR(NSS,NSS),USOI(NSS,NSS)
      DIMENSION DENSOUT(6,NBTRI)
      Dimension IZMR(3),IZMI(3)
      Dimension IZMR2(3),IZMI2(3)
      DIMENSION IOFF(8)
      CHARACTER*8 CHARTYPE
      INTEGER ASS,BSS,ASF,BSF





C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

C The following creates an array that is used to
C map a specific spin state to the corresponding
C spin-free state and to its spin
C (see prprop.f and others)

      CALL GETMEM('MAPST','ALLO','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','ALLO','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','ALLO','INTE',LMAPMS,NSS)

      ISS=0
      DO ISF=1,NSTATE
        JOB=iWork(lJBNUM+ISF-1)
        MPLET=MLTPLT(JOB)

        DO MSPROJ=-MPLET+1,MPLET-1,2
          ISS=ISS+1
          IWORK(LMAPST-1+ISS)=ISF
          IWORK(LMAPSP-1+ISS)=MPLET
          IWORK(LMAPMS-1+ISS)=MSPROJ
        END DO
      END DO


c Allocate some arrays
c LSDMXR, etc      DM/TDM for this iteration
c LSDMXR2, etc     Accumulated DM/TDM
c LTMPR,I          Temporary array for U*AU multiplication
c LTDMZZ           DM/TDM as read from file
c LSCR             Scratch for expansion of LTDMZZ
      CALL GETMEM('TSDMXR','ALLO','REAL',LSDMXR,NBTRI)
      CALL GETMEM('TSDMYR','ALLO','REAL',LSDMYR,NBTRI)
      CALL GETMEM('TSDMZR','ALLO','REAL',LSDMZR,NBTRI)
      CALL GETMEM('TSDMXI','ALLO','REAL',LSDMXI,NBTRI)
      CALL GETMEM('TSDMYI','ALLO','REAL',LSDMYI,NBTRI)
      CALL GETMEM('TSDMZI','ALLO','REAL',LSDMZI,NBTRI)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMXR),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMYR),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMZR),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMXI),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMYI),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMZI),1)
      IZMR(1)=LSDMXR
      IZMR(2)=LSDMYR
      IZMR(3)=LSDMZR
      IZMI(1)=LSDMXI
      IZMI(2)=LSDMYI
      IZMI(3)=LSDMZI

      CALL GETMEM('TSDMXR2','ALLO','REAL',LSDMXR2,NBTRI)
      CALL GETMEM('TSDMYR2','ALLO','REAL',LSDMYR2,NBTRI)
      CALL GETMEM('TSDMZR2','ALLO','REAL',LSDMZR2,NBTRI)
      CALL GETMEM('TSDMXI2','ALLO','REAL',LSDMXI2,NBTRI)
      CALL GETMEM('TSDMYI2','ALLO','REAL',LSDMYI2,NBTRI)
      CALL GETMEM('TSDMZI2','ALLO','REAL',LSDMZI2,NBTRI)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMXR2),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMYR2),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMZR2),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMXI2),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMYI2),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMZI2),1)
      IZMR2(1)=LSDMXR2
      IZMR2(2)=LSDMYR2
      IZMR2(3)=LSDMZR2
      IZMI2(1)=LSDMXI2
      IZMI2(2)=LSDMYI2
      IZMI2(3)=LSDMZI2

      CALL GETMEM('TSDMTMPR','ALLO','REAL',LTMPR,NBTRI)
      CALL GETMEM('TSDMTMPI','ALLO','REAL',LTMPI,NBTRI)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LTMPR),1)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LTMPI),1)

      CALL GETMEM('TDMSCR','ALLO','REAL',LSCR,NBTRI)
c zeroed inside the loop
c      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSCR),1)

      CALL GETMEM('TDMZZ','ALLO','REAL',LTDMZZ,NTDMZZ)
      CALL DCOPY_(NTDMZZ,0.0D00,0,WORK(LTDMZZ),1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MAIN LOOP OVER KSF/LSF
C WRITTEN AS IN PRPROP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CORRESPONDING SPIN-FREE STATES OF THE
C REQUESTED SPIN STATES
      ASF=IWORK(LMAPST-1+ASS)
      BSF=IWORK(LMAPST-1+BSS)

      DO KSS=1,NSS
       KSF=IWORK(LMAPST-1+KSS)
       MPLETK=IWORK(LMAPSP-1+KSS)
       MSPROJK=IWORK(LMAPMS-1+KSS)

       DO LSS=1,NSS
        LSF=IWORK(LMAPST-1+LSS)
        MPLETL=IWORK(LMAPSP-1+LSS)
        MSPROJL=IWORK(LMAPMS-1+LSS)

        JOB1=iWork(lJBNUM+KSF-1)
        JOB2=iWork(lJBNUM+LSF-1)
        LSYM1=IRREP(JOB1)
        LSYM2=IRREP(JOB2)
        ISY12=MUL(LSYM1,LSYM2)

C SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS
        IOF=0
        DO ISY1=1,NSYM
          ISY2=MUL(ISY1,ISY12)
          IF(ISY1.LT.ISY2) GOTO 10
          IOFF(ISY1)=IOF
          IOFF(ISY2)=IOF
          NB1=NBASF(ISY1)
          NB2=NBASF(ISY2)
          NB12=NB1*NB2
          IF(ISY1.EQ.ISY2) NB12=(NB12+NB1)/2
          IOF=IOF+NB12
  10      CONTINUE
        END DO


c These are going to be zero, so head them off at the pass
        IF(ITYPE.LE.2
     &     .AND.(MPLETK.NE.MPLETL.OR.MSPROJK.NE.MSPROJL)) GOTO 2200


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C Transition density matrices, TDMZZ, in AO basis.
C WDMZZ similar, but WE-reduced 'triplet' densities.
C TDMZZ will store either, depending on the type
        CALL DCOPY_(NTDMZZ,0.0D00,0,WORK(LTDMZZ),1)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C IDTDM: TOC array for transition 1-matrices
c TDMZZ is stored on disk from i = 1, NSTATE j=1, i
c so swap if needed
        IF(LSF.GT.KSF) THEN
          IDISK=iWork(lIDTDM+(LSF-1)*NSTATE+KSF-1)
        ELSE
          IDISK=iWork(lIDTDM+(KSF-1)*NSTATE+LSF-1)
        END IF
        CALL DDAFILE(LUTDM,2,WORK(LTDMZZ),NTDMZZ,IDISK)

c I Don't know what is stored between TDMZZ and WDMZZ,
c but store it in TDMZZ then overwrite
c (see mectl.f)
        IF(ITYPE.GE.3) THEN
          CALL DCOPY_(NTDMZZ,0.0D00,0,WORK(LTDMZZ),1)
          CALL DDAFILE(LUTDM,2,WORK(LTDMZZ),NTDMZZ,IDISK)
          CALL DCOPY_(NTDMZZ,0.0D00,0,WORK(LTDMZZ),1)
          CALL DDAFILE(LUTDM,2,WORK(LTDMZZ),NTDMZZ,IDISK)
C NOTE-the TD matrix as read in has an incorrect sign
          CALL DSCAL_(NTDMZZ,-1.0d0,WORK(LTDMZZ),1)
        END IF


c Anti-hermitian properties need a little fixing
        IF((ITYPE.EQ.2.OR.ITYPE.EQ.4).AND.(KSF.LE.LSF))
     &          CALL DSCAL_(NTDMZZ,-1.0d0,WORK(LTDMZZ),1)


C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
        CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSCR),1)

C This code expands the NTDMZZ-size matrix into
C an NBTRI-sized matrix
C SPECIAL CASE: DIAGONAL SYMMETRY BLOCKS.
C NOTE: During the folding, all off-diagonal values
C  are doubled. This factor of 2 is made up when
C  multiplying just the triangle with the AO matrix (DDOT)
C  rather than the entire matrix

C TODO - get rid of if statements in inner loop
ccccccccccccccccccc

c DIAGONAL SYMMETRY BLOCKS
        IF(ISY12.EQ.1) THEN
          IOF=0
          ITD=0
          DO 100 ISY=1,NSYM
            NB=NBASF(ISY)
            IF(NB.EQ.0) GOTO 100
            DO 90 J=1,NB
              DO 90 I=1,NB
                ITD=ITD+1
                TDM=WORK(LTDMZZ-1+ITD)
                IF(I.GE.J) THEN
                  IJ=IOF+(I*(I-1))/2+J
                  IF(I.GT.J) THEN
                    IF(ITYPE.EQ.2) WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)+TDM
                    IF(ITYPE.EQ.4) WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)+TDM
                  END IF
                ELSE
                  IJ=IOF+(J*(J-1))/2+I
                  IF(ITYPE.EQ.2) WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)-TDM
                  IF(ITYPE.EQ.4) WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)-TDM
                END IF
                IF(ITYPE.EQ.1) WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)+TDM
                IF(ITYPE.EQ.3) WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)+TDM
90          CONTINUE
            IOF=IOF+(NB*(NB+1))/2
100       CONTINUE
        ELSE
C GENERAL CASE, NON-DIAGONAL SYMMETRY BLOCKS
C THEN LOOP OVER ELEMENTS OF TDMZZ
          ITD=0
          DO 200 ISY1=1,NSYM
            NB1=NBASF(ISY1)
            IF(NB1.EQ.0) GOTO 200
            ISY2=MUL(ISY1,ISY12)
            NB2=NBASF(ISY2)
            IF(NB2.EQ.0) GOTO 200
            IF(ISY1.GT.ISY2) THEN
              DO 180 J=1,NB2
                DO 180 I=1,NB1
                  ITD=ITD+1
                  TDM=WORK(LTDMZZ-1+ITD)
                  IJ=IOFF(ISY1)+I+NB1*(J-1)
                  WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)+TDM
180           CONTINUE
            ELSE
              DO 190 J=1,NB2
                DO 190 I=1,NB1
                  ITD=ITD+1
                  TDM=WORK(LTDMZZ-1+ITD)
                  IJ=IOFF(ISY2)+J+NB2*(I-1)
                  WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)-TDM
190           CONTINUE
            END IF
200       CONTINUE
        END IF


c ie, see how AMFI is processed in soeig.f
        CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMXR),1)
        CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMYR),1)
        CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMZR),1)
        CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMXI),1)
        CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMYI),1)
        CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSDMZI),1)

        IF(ITYPE.GE.3) THEN
          S1=0.5D0*DBLE(MPLETK-1)
          SM1=0.5D0*DBLE(MSPROJK)
          S2=0.5D0*DBLE(MPLETL-1)
          SM2=0.5D0*DBLE(MSPROJL)
          FACT=1.0D0/SQRT(DBLE(MPLETK))
          IF(MPLETK.EQ.MPLETL-2) FACT=-FACT

          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)
        END IF

        IF((ITYPE.EQ.1.OR.ITYPE.EQ.2)
     &          .AND.MPLETK.EQ.MPLETL
     &          .AND.MSPROJK.EQ.MSPROJL) THEN
          CALL DAXPY_(NBTRI,1.0d0,WORK(LSCR),1,WORK(LSDMZR),1)
        ELSE IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
          CALL DAXPY_(NBTRI,CGX,WORK(LSCR),1,WORK(LSDMXR),1)
          CALL DAXPY_(NBTRI,CGY,WORK(LSCR),1,WORK(LSDMYI),1)
          CALL DAXPY_(NBTRI,CG0,WORK(LSCR),1,WORK(LSDMZR),1)
        END IF


cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c SPINORBIT
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Sign of the left-hand imaginary part is handled
c when doing DAXPY
        URR=USOR(LSS,BSS)
        UIR=USOI(LSS,BSS)
        URL=USOR(KSS,ASS)
        UIL=USOI(KSS,ASS)

        DO IDIR=1,3
          CALL DCOPY_(NBTRI,0.0D00,0,WORK(LTMPR),1)
          CALL DCOPY_(NBTRI,0.0D00,0,WORK(LTMPI),1)

C right side
          CALL DAXPY_(NBTRI,       URR,WORK(IZMR(IDIR)),1,WORK(LTMPR),1)
          CALL DAXPY_(NBTRI,-1.0d0*UIR,WORK(IZMI(IDIR)),1,WORK(LTMPR),1)
          CALL DAXPY_(NBTRI,       UIR,WORK(IZMR(IDIR)),1,WORK(LTMPI),1)
          CALL DAXPY_(NBTRI,       URR,WORK(IZMI(IDIR)),1,WORK(LTMPI),1)

C left side
         CALL DAXPY_(NBTRI,       URL,WORK(LTMPR),1,WORK(IZMR2(IDIR)),1)
         CALL DAXPY_(NBTRI,       UIL,WORK(LTMPI),1,WORK(IZMR2(IDIR)),1)
         CALL DAXPY_(NBTRI,       URL,WORK(LTMPI),1,WORK(IZMI2(IDIR)),1)
         CALL DAXPY_(NBTRI,-1.0d0*UIL,WORK(LTMPR),1,WORK(IZMI2(IDIR)),1)
        END DO
cccccccccccccccccccccc
c END SPINORBIT STUFF
cccccccccccccccccccccc



C END MAIN LOOP OVER STATES (KSS,LSS)
 2200   CONTINUE

       END DO
      END DO


C Store this density to DENSOUT
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
        DO I=1,NBTRI
          DENSOUT(1,I)=WORK(LSDMXR2-1+I)
          DENSOUT(2,I)=WORK(LSDMYR2-1+I)
          DENSOUT(3,I)=WORK(LSDMZR2-1+I)
          DENSOUT(4,I)=WORK(LSDMXI2-1+I)
          DENSOUT(5,I)=WORK(LSDMYI2-1+I)
          DENSOUT(6,I)=WORK(LSDMZI2-1+I)
        END DO
      ELSE
        DO I=1,NBTRI
          DENSOUT(1,I)=WORK(LSDMZR2-1+I)
          DENSOUT(2,I)=WORK(LSDMZR2-1+I)
          DENSOUT(3,I)=WORK(LSDMZR2-1+I)
          DENSOUT(4,I)=WORK(LSDMZI2-1+I)
          DENSOUT(5,I)=WORK(LSDMZI2-1+I)
          DENSOUT(6,I)=WORK(LSDMZI2-1+I)
        END DO
      END IF

c Free memory
      CALL GETMEM('TDMSCR','FREE','REAL',LSCR,NBTRI)
      CALL GETMEM('TDMZZ','FREE','REAL',LTDMZZ,NTDMZZ)

      CALL GETMEM('TSDMXR','FREE','REAL',LSDMXR,NBTRI)
      CALL GETMEM('TSDMYR','FREE','REAL',LSDMYR,NBTRI)
      CALL GETMEM('TSDMZR','FREE','REAL',LSDMZR,NBTRI)
      CALL GETMEM('TSDMXI','FREE','REAL',LSDMXI,NBTRI)
      CALL GETMEM('TSDMYI','FREE','REAL',LSDMYI,NBTRI)
      CALL GETMEM('TSDMZI','FREE','REAL',LSDMZI,NBTRI)

      CALL GETMEM('TSDMTMPR','FREE','REAL',LTMPR,NBTRI)
      CALL GETMEM('TSDMTMPI','FREE','REAL',LTMPI,NBTRI)

      CALL GETMEM('TSDMXR2','FREE','REAL',LSDMXR2,NBTRI)
      CALL GETMEM('TSDMYR2','FREE','REAL',LSDMYR2,NBTRI)
      CALL GETMEM('TSDMZR2','FREE','REAL',LSDMZR2,NBTRI)
      CALL GETMEM('TSDMXI2','FREE','REAL',LSDMXI2,NBTRI)
      CALL GETMEM('TSDMYI2','FREE','REAL',LSDMYI2,NBTRI)
      CALL GETMEM('TSDMZI2','FREE','REAL',LSDMZI2,NBTRI)

      CALL GETMEM('MAPST','FREE','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','FREE','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','FREE','INTE',LMAPMS,NSS)

      RETURN
      END







      SUBROUTINE SONATORB_INT(DENS, CHARPROP, IC, CHARTYPE,ASS,BSS,NSS,
     &                        PROPVALXR,PROPVALYR,PROPVALZR,
     &                        PROPVALXI,PROPVALYI,PROPVALZI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SONATORB_INT')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION DENS(6,NBTRI)
      INTEGER ASS,BSS
      CHARACTER*8 CHARPROP, CHARTYPE
c      DIMENSION IOFF(8)




C SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS
c      IOF=0
c      DO ISY1=1,NSYM
c        ISY2=MUL(ISY1,ISY12)
c        IF(ISY1.LT.ISY2) GOTO 10
c        IOFF(ISY1)=IOF
c        IOFF(ISY2)=IOF
c        NB1=NBASF(ISY1)
c        NB2=NBASF(ISY2)
c        NB12=NB1*NB2
c        IF(ISY1.EQ.ISY2) NB12=(NB12+NB1)/2
c        IOF=IOF+NB12
c10      CONTINUE
c      END DO

C NOW DO INTEGRATION WITH AO MATRICES
C FOR THE EXPECTATION VALUE

C The following creates an array that is used to
C map a specific spin state to the corresponding
C spin-free state and to its spin
C (see prprop.f and others)

      CALL GETMEM('MAPST','ALLO','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','ALLO','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','ALLO','INTE',LMAPMS,NSS)

      ISS=0
      DO ISF=1,NSTATE
        JOB=iWork(lJBNUM+ISF-1)
        MPLET=MLTPLT(JOB)

        DO MSPROJ=-MPLET+1,MPLET-1,2
          ISS=ISS+1
          IWORK(LMAPST-1+ISS)=ISF
          IWORK(LMAPSP-1+ISS)=MPLET
          IWORK(LMAPMS-1+ISS)=MSPROJ
        END DO
      END DO


C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

C ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
c The extra 4 elements correspond to the nuclear contribution
c and the origin of the operator
      NIP=4+NBTRI
      CALL GETMEM('IP    ','ALLO','REAL',LIP,NIP)

C Get info from the stored integrals
c IOPT controls what is read.
c IOPT=1 Read the size information
c IOPT=0 Read the property
c IOPT=6 Read the property, skipping the nuclear contribution and the origin
c (see misc_util/OneFlags.fh)
      IOPT=1
      CALL iRDONE(IRC,IOPT,CHARPROP,IC,NSIZ,ISCHK)


C sanity check?
c      DO ISS=1,NSS
c        DO JSS=1,NSS
c          ISF=IWORK(LMAPST-1+ISS)
c          JSF=IWORK(LMAPST-1+JSS)
C COMBINED SYMMETRY OF STATES:
c          JOB1=JBNUM(ISF)
c          JOB2=JBNUM(JSF)
c          LSYM1=IRREP(JOB1)
c          LSYM2=IRREP(JOB2)
c          ISY12=MUL(LSYM1,LSYM2)
C THE SYMMETRY CHECK MASK:
c          MASK=2**(ISY12-1)
c          IF(MOD(ISCHK/MASK,2).EQ.0) GOTO ???
c        END DO
c      END DO

c Actually read the integral
      IOPT=0
      CALL RDONE(IRC,IOPT,CHARPROP,IC,WORK(LIP),ISCHK)

      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,'(6X,A)')'*** ERROR IN SUBROUTINE SONATORB ***'
        WRITE(6,'(6X,A)')'  FAILED IN READING FROM  ONEINT'
        WRITE(6,'(6X,A,A)')'  LABEL     = ',CHARPROP
        WRITE(6,'(6X,A,I2)')'  COMPONENT = ',IC
        WRITE(6,*)
        CALL ABEND()
      END IF

      PROPVALXR=0.0d0
      PROPVALYR=0.0d0
      PROPVALZR=0.0d0
      PROPVALXI=0.0d0
      PROPVALYI=0.0d0
      PROPVALZI=0.0d0

C The integral is NBTRI matrix
C The property is NBTRI matrix
c We only work with half the matrix. Therefore, this would
c have a factor of 2. However, the factor of 1/2 was missing
c in SONATORB.F from the symmetric/antisymmetric equations
      IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
        DO I=1,NBTRI
          PROPVALXR=PROPVALXR+WORK(LIP-1+I)*DENS(1,I)
          PROPVALYR=PROPVALYR+WORK(LIP-1+I)*DENS(2,I)
          PROPVALZR=PROPVALZR+WORK(LIP-1+I)*DENS(3,I)

          PROPVALXI=PROPVALXI+WORK(LIP-1+I)*DENS(4,I)
          PROPVALYI=PROPVALYI+WORK(LIP-1+I)*DENS(5,I)
          PROPVALZI=PROPVALZI+WORK(LIP-1+I)*DENS(6,I)
        END DO
      ELSE
        DO I=1,NBTRI
          PROPVALXI=PROPVALXI+WORK(LIP-1+I)*DENS(1,I)
          PROPVALYI=PROPVALYI+WORK(LIP-1+I)*DENS(2,I)
          PROPVALZI=PROPVALZI+WORK(LIP-1+I)*DENS(3,I)

          PROPVALXR=PROPVALXR-WORK(LIP-1+I)*DENS(4,I)
          PROPVALYR=PROPVALYR-WORK(LIP-1+I)*DENS(5,I)
          PROPVALZR=PROPVALZR-WORK(LIP-1+I)*DENS(6,I)
        END DO
      END IF


      IF(IPGLOB.GE.VERBOSE) THEN
      IF(ITYPE.EQ.1.OR.ITYPE.EQ.2) THEN
        WRITE(6,*)
        WRITE(6,*) "************************************"
        WRITE(6,*) "SONATORB EXPECTATION VALUES"
        WRITE(6,*) " PROPERTY: ", CHARPROP
        WRITE(6,*) " COMPONENT: ", IC
        WRITE(6,*) " TYPE: ", CHARTYPE
        WRITE(6,*) " STATE (K,L): ",ASS,BSS
        WRITE(6,*) "************************************"
        WRITE(6,*) "Property: Real: ",PROPVALZR
        WRITE(6,*) "Property: Imag: ",PROPVALZI
        WRITE(6,*) "************************************"
      ELSE
        WRITE(6,*)
        WRITE(6,*) "************************************"
        WRITE(6,*) "SONATORB EXPECTATION VALUES"
        WRITE(6,*) " PROPERTY: ", CHARPROP
        WRITE(6,*) " COMPONENT: ", IC
        WRITE(6,*) " TYPE: ", CHARTYPE
        WRITE(6,*) " STATE (K,L): ",ASS,BSS
        WRITE(6,*) "************************************"
        WRITE(6,*) "Property: Re(X): ",PROPVALXR
        WRITE(6,*) "Property: Re(Y): ",PROPVALYR
        WRITE(6,*) "Property: Re(Z): ",PROPVALZR
        WRITE(6,*) "Property: Im(X): ",PROPVALXI
        WRITE(6,*) "Property: Im(Y): ",PROPVALYI
        WRITE(6,*) "Property: Im(Z): ",PROPVALZI
        WRITE(6,*) "************************************"
      END IF
      END IF

c Free up un-needed space
      CALL GETMEM('IP    ','FREE','REAL',LIP,NIP)
      CALL GETMEM('MAPST','FREE','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','FREE','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','FREE','INTE',LMAPMS,NSS)

      RETURN
      END





      SUBROUTINE SONATORB_PLOT (DENS, FILEBASE, CHARTYPE, ASS, BSS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SONATORB_PLOT')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION DENS(6,NBTRI)
      CHARACTER*25 FNAME
      CHARACTER(LEN=*) FILEBASE
      CHARACTER*16 KNUM
      CHARACTER*16 FNUM,XNUM
      CHARACTER*8 CHARTYPE
      CHARACTER CDIR
      INTEGER ASS,BSS




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTTING SECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c LSZZ  - AO Overlap integral
c LVEC  - AO Overlap eigenvectors
c LEIG  - AO Overlap eigenvalues
c LVEC2 - Eigenvectors of density matrix
c LSCR  - Temporary for matrix multiplication
C NOTE: LSCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL GETMEM('SZZ   ','ALLO','REAL',LSZZ,NBTRI)
      CALL GETMEM('VEC   ','ALLO','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','ALLO','REAL',LVEC2,NBMX2)
      CALL GETMEM('SCR   ','ALLO','REAL',LSCR,NBMX2)
      CALL GETMEM('EIG   ','ALLO','REAL',LEIG,NBST)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSZZ),1)
      CALL DCOPY_(NBSQ,0.0D00,0,WORK(LVEC),1)
      CALL DCOPY_(NBMX2,0.0D00,0,WORK(LVEC2),1)
      CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCR),1)
      CALL DCOPY_(NBST,0.0D00,0,WORK(LEIG),1)

      CALL GETMEM('VNAT  ','ALLO','REAL',LVNAT,NBSQ)
      CALL GETMEM('OCC   ','ALLO','REAL',LOCC,NBST)
      CALL DCOPY_(NBSQ,0.0D00,0,WORK(LVNAT),1)
      CALL DCOPY_(NBST,0.0D00,0,WORK(LOCC),1)

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=6
      ICMP=1
      ISYLAB=1
      CALL RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSZZ),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=LSZZ
      LV=LVEC
      LE=LEIG
      CALL VCLR(WORK(LVEC),1,NBSQ)
      DO 700 ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO 620 I=1,NB**2,(NB+1)
          WORK(LV-1+I)=1.0D00
620      CONTINUE
        CALL JACOB(WORK(LS),WORK(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO 630 I=1,NB
          EIG=WORK(LS1)
          WORK(LE1)=EIG
          X=1.0D00/SQRT(MAX(EIG,1.0D-14))
          CALL DSCAL_(NB,X,WORK(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
630      CONTINUE
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
700   CONTINUE

      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NBTRI)

      CALL GETMEM('TDMAT ','ALLO','REAL',LDMAT,NBMX2)

      IF(ITYPE.LE.2) THEN
        ISTART=3
        IEND=3
      ELSE
        ISTART=1
        IEND=3
      END IF

      DO IDIR=ISTART,IEND

        CDIR='?'
        IF(IDIR.EQ.1) CDIR='X'
        IF(IDIR.EQ.2) CDIR='Y'
        IF(IDIR.EQ.3) CDIR='Z'

        INV=1
        II2=0
        IOCC=0
        LV=LVEC
        LE=LEIG
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) GOTO 1750

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          CALL DCOPY_(NBMX2,0.0D0,0,WORK(LDMAT),1)
          CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCR),1)
          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)/2.0d0
              WORK(LDMAT-1+JI)=DENS(IDIR,II2)/2.0d0
            ELSE
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)
              WORK(LDMAT-1+JI)=DENS(IDIR,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LDMAT),NB,WORK(LV),NB,
     &                 0.0D0,WORK(LSCR),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LSCR),NB,
     &                 0.0D0,WORK(LDMAT),NB)

          ID1=1
          ID2=1
          DO I=1,NB
            EIG=WORK(LE-1+I)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID1),NB)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCR),1)
          ISCR=LSCR
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              WORK(ISCR)=(WORK(LDMAT-1+IJ)+WORK(LDMAT-1+JI))/2.0d0

c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) WORK(ISCR)=WORK(ISCR)*2.0d0
              ISCR=ISCR+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          CALL DCOPY_(NBMX2,0.0D0,0,WORK(LVEC2),1)
          CALL DCOPY_(NB,1.0D0,0,WORK(LVEC2),NB+1)

          CALL JACOB(WORK(LSCR),WORK(LVEC2),NB,NB)
          CALL JACORD(WORK(LSCR),WORK(LVEC2),NB,NB)

C JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=LSCR-1
          DO I=1,NB
            II=II+I
            WORK(LOCC-1+IOCC+NB+1-I)=WORK(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVALUES IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LVEC2),NB,
     &                 0.0D0,WORK(LSCR),NB)
          I1=LSCR
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,WORK(I1),1,WORK(LVNAT-1+I2),1)
            I1=I1+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
1750      CONTINUE
        END DO

C WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2)') ASS,".",BSS
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".",CDIR
       END IF
       WRITE(FNUM,'(I8)') BSS
       FNUM=ADJUSTL(FNUM)
       IF (ASS.NE.BSS) THEN
         WRITE(XNUM,'(I8,A)') ASS,'_'//TRIM(FNUM)
         FNUM=ADJUSTL(XNUM)
       END IF
       IF (ITYPE.GT.2) FNUM=CDIR//TRIM(FNUM)

       FNAME=FILEBASE//'.'//TRIM(FNUM)
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        IFOCC=1
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &     WORK(LVNAT), WORK(LOCC), Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_PLOT", WORK(LVNAT), 1, 4)

c    ONLYFOR NATURAL ORBITALS
      if(ITYPE.EQ.1)
     &       CALL ADD_INFO("SONATORB_NO_OCC", WORK(LOCC), NBASF, 4)

      END DO

      CALL GETMEM('TDMAT ','FREE','REAL',LDMAT,NBMX2)
      CALL GETMEM('VEC   ','FREE','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','FREE','REAL',LVEC2,NBMX2)
      CALL GETMEM('SCR   ','FREE','REAL',LSCR,NBMX2)
      CALL GETMEM('EIG   ','FREE','REAL',LEIG,NBST)
      CALL GETMEM('VNAT  ','FREE','REAL',LVNAT,NBSQ)
      CALL GETMEM('OCC   ','FREE','REAL',LOCC,NBST)


      RETURN
      END




      SUBROUTINE SONATORB_CPLOT (DENS, FILEBASE, CHARTYPE, ASS, BSS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SONATORB_CPLOT')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION DENS(6,NBTRI)
      CHARACTER*25 FNAME
      CHARACTER(LEN=*) FILEBASE
      CHARACTER*16 KNUM
      CHARACTER*16 FNUM,XNUM
      CHARACTER*8 CHARTYPE
      CHARACTER CDIR
      INTEGER ASS,BSS




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTTING SECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c LSZZ  - AO Overlap integral
c LVEC  - AO Overlap eigenvectors
c LEIG  - AO Overlap eigenvalues
c LVEC2 - Eigenvectors of density matrix
c LSCR  - Temporary for matrix multiplication
C NOTE: LSCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL GETMEM('SZZ   ','ALLO','REAL',LSZZ,NBTRI)
      CALL GETMEM('VEC   ','ALLO','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','ALLO','REAL',LVEC2,NBMX2)
      CALL GETMEM('VEC2I  ','ALLO','REAL',LVEC2I,NBMX2)
      CALL GETMEM('SCR   ','ALLO','REAL',LSCR,NBMX2)
      CALL GETMEM('SCRI   ','ALLO','REAL',LSCRI,NBMX2)
      CALL GETMEM('EIG   ','ALLO','REAL',LEIG,NBST)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSZZ),1)
      CALL DCOPY_(NBSQ,0.0D00,0,WORK(LVEC),1)
      CALL DCOPY_(NBMX2,0.0D00,0,WORK(LVEC2),1)
      CALL DCOPY_(NBMX2,0.0D00,0,WORK(LVEC2I),1)
      CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCR),1)
      CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCRI),1)
      CALL DCOPY_(NBST,0.0D00,0,WORK(LEIG),1)

      CALL GETMEM('VNAT  ','ALLO','REAL',LVNAT,NBSQ)
      CALL GETMEM('VNATI  ','ALLO','REAL',LVNATI,NBSQ)
      CALL GETMEM('OCC   ','ALLO','REAL',LOCC,NBST)
      CALL DCOPY_(NBSQ,0.0D00,0,WORK(LVNAT),1)
      CALL DCOPY_(NBSQ,0.0D00,0,WORK(LVNATI),1)
      CALL DCOPY_(NBST,0.0D00,0,WORK(LOCC),1)

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=6
      ICMP=1
      ISYLAB=1
      CALL RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSZZ),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=LSZZ
      LV=LVEC
      LE=LEIG
      CALL VCLR(WORK(LVEC),1,NBSQ)
      DO 1700 ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO 1620 I=1,NB**2,(NB+1)
          WORK(LV-1+I)=1.0D00
1620      CONTINUE
        CALL JACOB(WORK(LS),WORK(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO 1630 I=1,NB
          EIG=WORK(LS1)
          WORK(LE1)=EIG
          X=1.0D00/SQRT(MAX(EIG,1.0D-14))
          CALL DSCAL_(NB,X,WORK(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
1630      CONTINUE
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
1700   CONTINUE

      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NBTRI)

      CALL GETMEM('TDMAT ','ALLO','REAL',LDMAT,NBMX2)
      CALL GETMEM('TDMATI ','ALLO','REAL',LDMATI,NBMX2)

      IF(ITYPE.LE.2) THEN
        ISTART=3
        IEND=3
      ELSE
        ISTART=1
        IEND=3
      END IF

      DO IDIR=ISTART,IEND

          CDIR='?'
          IF(IDIR.EQ.1) CDIR='X'
          IF(IDIR.EQ.2) CDIR='Y'
          IF(IDIR.EQ.3) CDIR='Z'


cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
C read in ao matrix for angmom or mltpl
      CALL GETMEM('SANG  ','ALLO','REAL',LSANG,NBTRI)
      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSANG),1)

      IRC=-1
      IOPT=6

      IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
        ICMP=1
        CALL iRDONE(IRC,   1,'MLTPL  0',ICMP,NSIZ,       ISYLAB)
        CALL  RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSANG),ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(6,*)
          WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(6,*)'      MLTPL0 INTEGRALS ARE NOT AVAILABLE'
          WRITE(6,*)'      IRC:',IRC
          WRITE(6,*)
          CALL ABEND()
        END IF

      ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
        ICMP=3
        CALL iRDONE(IRC,   1,'ANGMOM  ',ICMP,NSIZ,       ISYLAB)
        CALL  RDONE(IRC,IOPT,'ANGMOM  ',ICMP,WORK(LSANG),ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(6,*)
          WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(6,*)'      ANGMOM INTEGRALS ARE NOT AVAILABLE'
          WRITE(6,*)'      IRC:',IRC
          WRITE(6,*)
          CALL ABEND()
        END IF

      END IF

cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
        INV=1
        II2=0
        IOCC=0
        LV=LVEC
        LE=LEIG
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) GOTO 1800

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          CALL DCOPY_(NBMX2,0.0D0,0,WORK(LDMAT),1)
          CALL DCOPY_(NBMX2,0.0D0,0,WORK(LDMATI),1)
          CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCR),1)
          CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCRI),1)

          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)/2.0d0
              WORK(LDMAT-1+JI)=DENS(IDIR,II2)/2.0d0
              WORK(LDMATI-1+IJ)=-1.0d0*DENS(IDIR+3,II2)/2.0d0
              WORK(LDMATI-1+JI)=DENS(IDIR+3,II2)/2.0d0
            ELSE
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)
              WORK(LDMATI-1+JI)=DENS(IDIR+3,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LDMAT),NB,WORK(LV),NB,
     &                 0.0D0,WORK(LSCR),NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LDMATI),NB,WORK(LV),NB,
     &                 0.0D0,WORK(LSCRI),NB)



          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LSCR),NB,
     &                 0.0D0,WORK(LDMAT),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LSCRI),NB,
     &                 0.0D0,WORK(LDMATI),NB)

          ID1=1
          ID2=1
          DO I=1,NB
            EIG=WORK(LE-1+I)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID1),NB)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID2),1)
            CALL DSCAL_(NB,EIG,WORK(LDMATI-1+ID1),NB)
            CALL DSCAL_(NB,EIG,WORK(LDMATI-1+ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCR),1)
          CALL DCOPY_(NBMX2,0.0D00,0,WORK(LSCRI),1)

          ISCR=LSCR
          ISCRI=LSCRI
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              WORK(ISCR)=(WORK(LDMAT-1+JI)+WORK(LDMAT-1+IJ))/2.0d0
              WORK(ISCRI)=(WORK(LDMATI-1+JI)-WORK(LDMATI-1+IJ))/2.0d0
c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) WORK(ISCR)=WORK(ISCR)*2.0d0
              IF(ITYPE.GE.3) WORK(ISCRI)=WORK(ISCRI)*2.0d0
              ISCR=ISCR+1
              ISCRI=ISCRI+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          CALL DCOPY_(NBMX2,0.0D0,0,WORK(LVEC2),1)
          CALL DCOPY_(NBMX2,0.0D0,0,WORK(LVEC2I),1)

          CALL CPLOT_DIAG(WORK(LSCR),WORK(LSCRI), NB,
     &                    WORK(LVEC2),WORK(LVEC2I))

C LAPACK ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=LSCR-1
          DO I=1,NB
            II=II+I
            WORK(LOCC-1+IOCC+NB+1-I)=WORK(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LVEC2),NB,
     &                 0.0D0,WORK(LSCR),NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LVEC2I),NB,
     &                 0.0D0,WORK(LSCRI),NB)

          I1=LSCR
          I1I=LSCRI
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,WORK(I1),1,WORK(LVNAT-1+I2),1)
            CALL DCOPY_(NB,WORK(I1I),1,WORK(LVNATI-1+I2),1)
            I1=I1+NB
            I1I=I1I+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
1800      CONTINUE
        END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC TESTING
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(IPGLOB.GE.DEBUG) THEN

      CALL GETMEM('SANGF ','ALLO','REAL',LSANGF,NBMX**2)
      CALL GETMEM('SANGTR  ','ALLO','REAL',LSANGTR,NBMX**2)
      CALL GETMEM('SANGTI  ','ALLO','REAL',LSANGTI,NBMX**2)
      CALL GETMEM('SANGTR2  ','ALLO','REAL',LSANGTR2,NBMX**2)
      CALL GETMEM('SANGTI2  ','ALLO','REAL',LSANGTI2,NBMX**2)
      CALL DCOPY_(NBMX**2,0.0D00,0,WORK(LSANGF),1)
      CALL DCOPY_(NBMX**2,0.0D00,0,WORK(LSANGTR),1)
      CALL DCOPY_(NBMX**2,0.0D00,0,WORK(LSANGTI),1)
      CALL DCOPY_(NBMX**2,0.0D00,0,WORK(LSANGTR2),1)
      CALL DCOPY_(NBMX**2,0.0D00,0,WORK(LSANGTI2),1)

      INV=0
      INV2=0
      II=0
      SUM = 0.0d0
      SUMI = 0.0d0

      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        IF(NB.EQ.0) GOTO 1860

c       Expand integrals for this symmetry to full storage
        CALL DCOPY_(NBMX**2,0.0d0,0,WORK(LSANGF),1)

        DO J=1,NB
        DO I=1,J
          IJ=NB*(J-1)+I-1
          JI=NB*(I-1)+J-1

          WORK(LSANGF+JI) = WORK(LSANG+II)

          IF(I.NE.J) THEN
            IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
              WORK(LSANGF+IJ) = 1.0d0 * WORK(LSANG+II)
            ELSE
              WORK(LSANGF+IJ) = WORK(LSANG+II)
            END IF
          END IF

          II=II+1

        END DO
        END DO

        IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LSANGF),NB,
     &             WORK(LVNAT+INV),NB,0.0d0,WORK(LSANGTR),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LSANGF),NB,
     &              WORK(LVNATI+INV),NB,0.0d0,WORK(LSANGTI),NB)

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTR2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTR2),NB)

          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTI2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTI2),NB)

        ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LSANGF),NB,
     &             WORK(LVNAT+INV),NB,0.0d0,WORK(LSANGTI),NB)
          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,WORK(LSANGF),NB,
     &             WORK(LVNATI+INV),NB,0.0d0,WORK(LSANGTR),NB)

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTR2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTR2),NB)

          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTI2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTI2),NB)

        END IF

c Sum over the trace
        DO I = 1,NB
          IJ = I+(I-1)*NB-1
          SUM  = SUM  + WORK(LOCC-1+I+INV2) * WORK(LSANGTR2+IJ)
          SUMI = SUMI + WORK(LOCC-1+I+INV2) * WORK(LSANGTI2+IJ)
        END DO

1860    CONTINUE

        INV=INV+NB**2
        INV2=INV2+NB

      END DO

        WRITE(6,*) "Ben P TEST for JA:"
        WRITE(6,*) "REAL: ",SUM
        WRITE(6,*) "IMAG: ",SUMI

        CALL GETMEM('SANGF ','FREE','REAL',LSANGF,NBMX**2)
        CALL GETMEM('SANGTR  ','FREE','REAL',LSANGTR,NBMX**2)
        CALL GETMEM('SANGTI  ','FREE','REAL',LSANGTI,NBMX**2)
        CALL GETMEM('SANGTR2  ','FREE','REAL',LSANGTR2,NBMX**2)
        CALL GETMEM('SANGTI2  ','FREE','REAL',LSANGTI2,NBMX**2)
      END IF ! IPGLOB >= DEBUG

      CALL GETMEM('SANG  ','FREE','REAL',LSANG,NBTRI)

C WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
C REAL PART
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".","R"
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A,A,A)')ASS,".",BSS,".",CDIR,".","R"
       END IF
       WRITE(FNUM,'(I8)') BSS
       FNUM=ADJUSTL(FNUM)
       IF (ASS.NE.BSS) THEN
         WRITE(XNUM,'(I8,A)') ASS,'_'//TRIM(FNUM)
         FNUM=ADJUSTL(XNUM)
       END IF
       IF (ITYPE.GT.2) FNUM=CDIR//TRIM(FNUM)

       FNAME=FILEBASE//'.'//TRIM(FNUM)//'.R'
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        IFOCC=1
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &     WORK(LVNAT), WORK(LOCC), Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

C IMAGINARY PART
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".","I"
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A,A,A)')ASS,".",BSS,".",CDIR,".","I"
       END IF

       FNAME=FILEBASE//'.'//TRIM(FNUM)//'.I'
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        IFOCC=1
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &     WORK(LVNATI), WORK(LOCC), Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_CPLOTR", WORK(LVNAT), 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTI", WORK(LVNATI), 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTO", WORK(LOCC), 1, 4)

      END DO

      CALL GETMEM('TDMAT ','FREE','REAL',LDMAT,NBMX2)
      CALL GETMEM('TDMATI ','FREE','REAL',LDMATI,NBMX2)
      CALL GETMEM('VEC   ','FREE','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','FREE','REAL',LVEC2,NBMX2)
      CALL GETMEM('VEC2I  ','FREE','REAL',LVEC2I,NBMX2)
      CALL GETMEM('SCR   ','FREE','REAL',LSCR,NBMX2)
      CALL GETMEM('SCRI   ','FREE','REAL',LSCRI,NBMX2)
      CALL GETMEM('EIG   ','FREE','REAL',LEIG,NBST)
      CALL GETMEM('VNAT  ','FREE','REAL',LVNAT,NBSQ)
      CALL GETMEM('VNATI  ','FREE','REAL',LVNATI,NBSQ)
      CALL GETMEM('OCC   ','FREE','REAL',LOCC,NBST)

      RETURN
      END





      SUBROUTINE CPLOT_DIAG(MATR, MATI, DIM, EIGVECR, EIGVECI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER DIM
      REAL*8 MATR(DIM*(DIM+1)/2),MATI(DIM*(DIM+1)/2)
      REAL*8 EIGVECR(DIM,DIM),EIGVECI(DIM,DIM)
      REAL*8 CEIGVAL(DIM)
      COMPLEX*16 MATFULL((DIM*(DIM+1)/2))
      COMPLEX*16 CEIGVEC(DIM,DIM)
      COMPLEX*16 ZWORK(2*DIM-1)
      REAL*8 RWORK(3*DIM-2)
      INTEGER INFO

      DO J=1,(DIM*(DIM+1)/2)
          MATFULL(J) = DCMPLX(MATR(J),MATI(J))
c          MATFULL(J) = DCMPLX(MATR(J),0.0d0)
      END DO


      call zhpev_('V','U',DIM,MATFULL,CEIGVAL,
     &           CEIGVEC,DIM,ZWORK,RWORK,INFO)


      IF(INFO.NE.0) THEN
          WRITE(6,*) "Error in diagonalization"
          WRITE(6,*) "INFO: ",INFO
          CALL ABEND()
      END IF

      DO I=1,DIM
      DO J=1,DIM
          EIGVECR(I,J) = REAL(CEIGVEC(I,J))
          EIGVECI(I,J) = AIMAG(CEIGVEC(I,J))
      END DO
      END DO

      CALL DCOPY_(DIM*(DIM+1)/2,0.0D00,0,MATR,1)
      CALL DCOPY_(DIM*(DIM+1)/2,0.0D00,0,MATI,1)

      DO J=1,DIM
         MATR((J*(J-1)/2)+J) = CEIGVAL(J)
      END DO

      RETURN
      END
