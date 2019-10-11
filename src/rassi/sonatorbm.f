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
      SUBROUTINE SONATORBM(CHARTYPE,
     &                     USOR,USOI,ASS,BSS,NSS,
     &                     iOpt,ROTMAT,DENSOUT)
      use rassi_aux, only : idisk_TDM
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
      Dimension ROTMAT(3,3)
      INTEGER ASS,BSS,ASF,BSF




c VV: dummy initialization
      CGY=-1
      CGX=-1
      CG0=-1
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
        iEmpty=iDisk_TDM(LSF,KSF,2)
        IDISK=iDisk_TDM(LSF,KSF,1)
        iOpt=2
        IF (ITYPE.GE.3) Then
           iGo=4
           CALL dens2file(Work(LTDMZZ),Work(LTDMZZ),Work(LTDMZZ),
     &                    nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo)
C NOTE-the TD matrix as read in has an incorrect sign
           CALL DSCAL_(NTDMZZ,-1.0d0,WORK(LTDMZZ),1)
        Else
           iGo=1
           CALL dens2file(Work(LTDMZZ),Work(LTDMZZ),Work(LTDMZZ),
     &                    nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo)
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
          If (iOpt.eq.1) Then
          CALL DAXPY_(NBTRI,CGX*ROTMAT(1,1),WORK(LSCR),1,WORK(LSDMXR),1)
          CALL DAXPY_(NBTRI,CGY*ROTMAT(2,1),WORK(LSCR),1,WORK(LSDMXI),1)
          CALL DAXPY_(NBTRI,CG0*ROTMAT(3,1),WORK(LSCR),1,WORK(LSDMXR),1)

          CALL DAXPY_(NBTRI,CGX*ROTMAT(1,2),WORK(LSCR),1,WORK(LSDMYR),1)
          CALL DAXPY_(NBTRI,CGY*ROTMAT(2,2),WORK(LSCR),1,WORK(LSDMYI),1)
          CALL DAXPY_(NBTRI,CG0*ROTMAT(3,2),WORK(LSCR),1,WORK(LSDMYR),1)

          CALL DAXPY_(NBTRI,CGX*ROTMAT(1,3),WORK(LSCR),1,WORK(LSDMZR),1)
          CALL DAXPY_(NBTRI,CGY*ROTMAT(2,3),WORK(LSCR),1,WORK(LSDMZI),1)
          CALL DAXPY_(NBTRI,CG0*ROTMAT(3,3),WORK(LSCR),1,WORK(LSDMZR),1)
          Else
          CALL DAXPY_(NBTRI,CGX,WORK(LSCR),1,WORK(LSDMXR),1)
          CALL DAXPY_(NBTRI,CGY,WORK(LSCR),1,WORK(LSDMYI),1)
          CALL DAXPY_(NBTRI,CG0,WORK(LSCR),1,WORK(LSDMZR),1)
          End If
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
