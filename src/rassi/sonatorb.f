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
      SUBROUTINE SONATORB(CHARTYPE,
     &                   USOR,USOI,ASS,BSS,NSS,
     &                   DENSOUT)
      use rassi_aux, only : iDisk_TDM
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
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMXR),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMYR),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMZR),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMXI),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMYI),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMZI),1)
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
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMXR2),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMYR2),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMZR2),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMXI2),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMYI2),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMZI2),1)
      IZMR2(1)=LSDMXR2
      IZMR2(2)=LSDMYR2
      IZMR2(3)=LSDMZR2
      IZMI2(1)=LSDMXI2
      IZMI2(2)=LSDMYI2
      IZMI2(3)=LSDMZI2

      CALL GETMEM('TSDMTMPR','ALLO','REAL',LTMPR,NBTRI)
      CALL GETMEM('TSDMTMPI','ALLO','REAL',LTMPI,NBTRI)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LTMPR),1)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LTMPI),1)

      CALL GETMEM('TDMSCR','ALLO','REAL',LSCR,NBTRI)
c zeroed inside the loop
c      CALL DCOPY_(NBTRI,0.0D00,0,WORK(LSCR),1)

      CALL GETMEM('TDMZZ','ALLO','REAL',LTDMZZ,NTDMZZ)
      CALL DCOPY_(NTDMZZ,[0.0D00],0,WORK(LTDMZZ),1)

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
        Call mk_IOFF(IOFF,nSYM,NBASF,ISY12)

c These are going to be zero, so head them off at the pass
        IF(ITYPE.LE.2
     &     .AND.(MPLETK.NE.MPLETL.OR.MSPROJK.NE.MSPROJL)) GOTO 2200


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C Transition density matrices, TDMZZ, in AO basis.
C WDMZZ similar, but WE-reduced 'triplet' densities.
C TDMZZ will store either, depending on the type
        CALL DCOPY_(NTDMZZ,[0.0D00],0,WORK(LTDMZZ),1)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C IDTDM: TOC array for transition 1-matrices
c TDMZZ is stored on disk from i = 1, NSTATE j=1, i
c so swap if needed
        IDISK=iDisk_TDM(KSF,LSF)
        CALL DDAFILE(LUTDM,2,WORK(LTDMZZ),NTDMZZ,IDISK)

c I Don't know what is stored between TDMZZ and WDMZZ,
c but store it in TDMZZ then overwrite
c (see mectl.f)
        IF(ITYPE.GE.3) THEN
          CALL DCOPY_(NTDMZZ,[0.0D00],0,WORK(LTDMZZ),1)
          CALL DDAFILE(LUTDM,2,WORK(LTDMZZ),NTDMZZ,IDISK)
          CALL DCOPY_(NTDMZZ,[0.0D00],0,WORK(LTDMZZ),1)
          CALL DDAFILE(LUTDM,2,WORK(LTDMZZ),NTDMZZ,IDISK)
C NOTE-the TD matrix as read in has an incorrect sign
          CALL DSCAL_(NTDMZZ,-1.0d0,WORK(LTDMZZ),1)
        END IF


c Anti-hermitian properties need a little fixing
        IF((ITYPE.EQ.2.OR.ITYPE.EQ.4).AND.(KSF.LE.LSF))
     &          CALL DSCAL_(NTDMZZ,-1.0d0,WORK(LTDMZZ),1)


C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
        CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSCR),1)

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
        CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMXR),1)
        CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMYR),1)
        CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMZR),1)
        CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMXI),1)
        CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMYI),1)
        CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSDMZI),1)

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
          CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LTMPR),1)
          CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LTMPI),1)

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
      DIMENSION IDUM(1)

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
      CALL iRDONE(IRC,IOPT,CHARPROP,IC,IDUM,ISCHK)
      IF (IRC.EQ.0) NSIZ=IDUM(1)

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
