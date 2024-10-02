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
      use rassi_global_arrays, only: JBNUM
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      CHARACTER(LEN=8) CHARTYPE
      INTEGER ASS,BSS,NSS
      Real*8 USOR(NSS,NSS),USOI(NSS,NSS)
      INTEGER iOpt
      Real*8 ROTMAT(3,3)
      Real*8 DENSOUT(6,NBTRI)

      Integer IZMR2(3),IZMI2(3)
      Integer IOFF(8)
      Integer, allocatable:: MAPST(:), MAPSP(:), MAPMS(:)
      Real*8, allocatable:: TMPR(:), TMPI(:), SCR(:), TDMZZ(:)
      Real*8, allocatable, Target:: SDMXR(:) ,SDMYR(:) ,SDMZR(:)
      Real*8, allocatable, Target:: SDMXI(:) ,SDMYI(:) ,SDMZI(:)
      Type A2_Array
        Real*8, Pointer:: A2(:)=>Null()
      End Type A2_Array
      Type (A2_array):: IZMR(3)
      Type (A2_array):: IZMI(3)

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

      CALL mma_allocate(MAPST,NSS,Label='MAPST')
      CALL mma_allocate(MAPSP,NSS,Label='MAPSP')
      CALL mma_allocate(MAPMS,NSS,Label='MAPMS')

      ISS=0
      DO ISF=1,NSTATE
        JOB=JBNUM(ISF)
        MPLET=MLTPLT(JOB)

        DO MSPROJ=-MPLET+1,MPLET-1,2
          ISS=ISS+1
          MAPST(ISS)=ISF
          MAPSP(ISS)=MPLET
          MAPMS(ISS)=MSPROJ
        END DO
      END DO


c Allocate some arrays
c SDMXR, etc      DM/TDM for this iteration
c SDMXR2, etc     Accumulated DM/TDM
c TMPR,I          Temporary array for U*AU multiplication
c TDMZZ           DM/TDM as read from file
c SCR             Scratch for expansion of TDMZZ
      CALL mma_allocate(SDMXR,NBTRI,Label='SDMXR')
      CALL mma_allocate(SDMYR,NBTRI,Label='SDMYR')
      CALL mma_allocate(SDMZR,NBTRI,Label='SDMZR')
      CALL mma_allocate(SDMXI,NBTRI,Label='SDMXI')
      CALL mma_allocate(SDMYI,NBTRI,Label='SDMYI')
      CALL mma_allocate(SDMZI,NBTRI,Label='SDMZI')
      SDMXR(:)=0.0D0
      SDMYR(:)=0.0D0
      SDMZR(:)=0.0D0
      SDMXI(:)=0.0D0
      SDMYI(:)=0.0D0
      SDMZI(:)=0.0D0
      IZMR(1)%A2=>SDMXR(:)
      IZMR(2)%A2=>SDMYR(:)
      IZMR(3)%A2=>SDMZR(:)
      IZMI(1)%A2=>SDMXI(:)
      IZMI(2)%A2=>SDMYI(:)
      IZMI(3)%A2=>SDMZI(:)

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

      CALL mma_allocate(TMPR,NBTRI,Label='TMPR')
      CALL mma_allocate(TMPI,NBTRI,Label='TMPI')
      TMPR(:)=0.0D0
      TMPI(:)=0.0D0

      CALL mma_allocate(SCR,NBTRI,Label='SCR')
c zeroed inside the loop
c      SCR(:)=0.0D0

      CALL mma_allocate(TDMZZ,NTDMZZ,Label='TDMZZ')
      TDMZZ(:)=0.0D0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MAIN LOOP OVER KSF/LSF
C WRITTEN AS IN PRPROP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CORRESPONDING SPIN-FREE STATES OF THE
C REQUESTED SPIN STATES
      DO KSS=1,NSS
       KSF=MAPST(KSS)
       MPLETK=MAPSP(KSS)
       MSPROJK=MAPMS(KSS)

       DO LSS=1,NSS
        LSF=MAPST(LSS)
        MPLETL=MAPSP(LSS)
        MSPROJL=MAPMS(LSS)

        JOB1=JBNUM(KSF)
        JOB2=JBNUM(LSF)
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
        CALL DCOPY_(NTDMZZ,[0.0D00],0,TDMZZ,1)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C IDTDM: TOC array for transition 1-matrices
c TDMZZ is stored on disk from i = 1, NSTATE j=1, i
c so swap if needed
        iEmpty=iDisk_TDM(KSF,LSF,2)
        IDISK=iDisk_TDM(KSF,LSF,1)
        iOpt=2
        IF (ITYPE.GE.3) Then
           iGo=4
           CALL dens2file(TDMZZ,TDMZZ,TDMZZ,
     &                    nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,KSF,LSF)
C NOTE-the TD matrix as read in has an incorrect sign
           CALL DSCAL_(NTDMZZ,-1.0d0,TDMZZ,1)
        Else
           iGo=1
           CALL dens2file(TDMZZ,TDMZZ,TDMZZ,
     &                    nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,KSF,LSF)
        END IF


c Anti-hermitian properties need a little fixing
        IF((ITYPE.EQ.2.OR.ITYPE.EQ.4).AND.(KSF.LE.LSF))
     &          CALL DSCAL_(NTDMZZ,-1.0d0,TDMZZ,1)


C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
        SCR(:)=0.0D0

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
          DO ISY=1,NSYM
            NB=NBASF(ISY)
            IF(NB.EQ.0) cycle
            DO J=1,NB
              DO I=1,NB
                ITD=ITD+1
                TDM=TDMZZ(ITD)
                IF(I.GE.J) THEN
                  IJ=IOF+(I*(I-1))/2+J
                  IF(I.GT.J) THEN
                    IF(ITYPE.EQ.2) SCR(IJ)=SCR(IJ)+TDM
                    IF(ITYPE.EQ.4) SCR(IJ)=SCR(IJ)+TDM
                  END IF
                ELSE
                  IJ=IOF+(J*(J-1))/2+I
                  IF(ITYPE.EQ.2) SCR(IJ)=SCR(IJ)-TDM
                  IF(ITYPE.EQ.4) SCR(IJ)=SCR(IJ)-TDM
                END IF
                IF(ITYPE.EQ.1) SCR(IJ)=SCR(IJ)+TDM
                IF(ITYPE.EQ.3) SCR(IJ)=SCR(IJ)+TDM
              END DO
            END DO
            IOF=IOF+(NB*(NB+1))/2
          END DO
        ELSE
C GENERAL CASE, NON-DIAGONAL SYMMETRY BLOCKS
C THEN LOOP OVER ELEMENTS OF TDMZZ
          ITD=0
          DO ISY1=1,NSYM
            NB1=NBASF(ISY1)
            IF(NB1.EQ.0) cycle
            ISY2=MUL(ISY1,ISY12)
            NB2=NBASF(ISY2)
            IF(NB2.EQ.0) cycle
            IF(ISY1.GT.ISY2) THEN
              DO J=1,NB2
                DO I=1,NB1
                  ITD=ITD+1
                  TDM=TDMZZ(ITD)
                  IJ=IOFF(ISY1)+I+NB1*(J-1)
                  SCR(IJ)=SCR(IJ)+TDM
                END DO
              END DO
            ELSE
              DO J=1,NB2
                DO I=1,NB1
                  ITD=ITD+1
                  TDM=TDMZZ(ITD)
                  IJ=IOFF(ISY2)+J+NB2*(I-1)
                  SCR(IJ)=SCR(IJ)-TDM
                END DO
              END DO
            END IF
          END DO
        END IF


c ie, see how AMFI is processed in soeig.f
        SDMXR(:)=0.0D0
        SDMYR(:)=0.0D0
        SDMZR(:)=0.0D0
        SDMXI(:)=0.0D0
        SDMYI(:)=0.0D0
        SDMZI(:)=0.0D0

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
          CALL DAXPY_(NBTRI,1.0d0,SCR,1,SDMZR,1)
        ELSE IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
          If (iOpt.eq.1) Then
          CALL DAXPY_(NBTRI,CGX*ROTMAT(1,1),SCR,1,SDMXR,1)
          CALL DAXPY_(NBTRI,CGY*ROTMAT(2,1),SCR,1,SDMXI,1)
          CALL DAXPY_(NBTRI,CG0*ROTMAT(3,1),SCR,1,SDMXR,1)

          CALL DAXPY_(NBTRI,CGX*ROTMAT(1,2),SCR,1,SDMYR,1)
          CALL DAXPY_(NBTRI,CGY*ROTMAT(2,2),SCR,1,SDMYI,1)
          CALL DAXPY_(NBTRI,CG0*ROTMAT(3,2),SCR,1,SDMYR,1)

          CALL DAXPY_(NBTRI,CGX*ROTMAT(1,3),SCR,1,SDMZR,1)
          CALL DAXPY_(NBTRI,CGY*ROTMAT(2,3),SCR,1,SDMZI,1)
          CALL DAXPY_(NBTRI,CG0*ROTMAT(3,3),SCR,1,SDMZR,1)
          Else
          CALL DAXPY_(NBTRI,CGX,SCR,1,SDMXR,1)
          CALL DAXPY_(NBTRI,CGY,SCR,1,SDMYI,1)
          CALL DAXPY_(NBTRI,CG0,SCR,1,SDMZR,1)
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
          CALL DCOPY_(NBTRI,[0.0D00],0,TMPR,1)
          CALL DCOPY_(NBTRI,[0.0D00],0,TMPI,1)

C right side
          CALL DAXPY_(NBTRI,       URR,IZMR(IDIR)%A2(:),1,TMPR,1)
          CALL DAXPY_(NBTRI,-1.0d0*UIR,IZMI(IDIR)%A2(:),1,TMPR,1)
          CALL DAXPY_(NBTRI,       UIR,IZMR(IDIR)%A2(:),1,TMPI,1)
          CALL DAXPY_(NBTRI,       URR,IZMI(IDIR)%A2(:),1,TMPI,1)

C left side
         CALL DAXPY_(NBTRI,       URL,TMPR,1,WORK(IZMR2(IDIR)),1)
         CALL DAXPY_(NBTRI,       UIL,TMPI,1,WORK(IZMR2(IDIR)),1)
         CALL DAXPY_(NBTRI,       URL,TMPI,1,WORK(IZMI2(IDIR)),1)
         CALL DAXPY_(NBTRI,-1.0d0*UIL,TMPR,1,WORK(IZMI2(IDIR)),1)
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
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(TDMZZ)

      Call mma_deallocate(SDMXR)
      Call mma_deallocate(SDMYR)
      Call mma_deallocate(SDMZR)
      Call mma_deallocate(SDMXI)
      Call mma_deallocate(SDMYI)
      Call mma_deallocate(SDMZI)
      IZMR(1)%A2=>Null()
      IZMR(2)%A2=>Null()
      IZMR(3)%A2=>Null()
      IZMI(1)%A2=>Null()
      IZMI(2)%A2=>Null()
      IZMI(3)%A2=>Null()

      Call mma_deallocate(TMPI)
      Call mma_deallocate(TMPR)

      CALL GETMEM('TSDMXR2','FREE','REAL',LSDMXR2,NBTRI)
      CALL GETMEM('TSDMYR2','FREE','REAL',LSDMYR2,NBTRI)
      CALL GETMEM('TSDMZR2','FREE','REAL',LSDMZR2,NBTRI)
      CALL GETMEM('TSDMXI2','FREE','REAL',LSDMXI2,NBTRI)
      CALL GETMEM('TSDMYI2','FREE','REAL',LSDMYI2,NBTRI)
      CALL GETMEM('TSDMZI2','FREE','REAL',LSDMZI2,NBTRI)

      Call mma_deallocate(MAPMS)
      Call mma_deallocate(MAPSP)
      Call mma_deallocate(MAPST)

      END SUBROUTINE SONATORBM
