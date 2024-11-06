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
* Copyright (C) 2021, Rulin Feng                                       *
************************************************************************
*       ***************************************************
*                Get TDM in AO basis from SO states
*       ****************************************************
*        This routine is modified from sonatorbm.f to give the
*        full transition density matrix (TDM) in atomic orbital
*        basis, which is not a nbtri sized matrix as in sonatorbm
*        but a full nbst**2 sized matrix.
*

      SUBROUTINE MAKETDMAO(CHARTYPE,
     &                     USOR,USOI,ASS,BSS,NSS,
     &                     iOpt,ROTMAT,DENSOUT)
      use rassi_aux, only : idisk_TDM
      use rassi_global_arrays, only: JBNUM
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, LSYM1, LSYM2, IRREP, MLTPLT
      use cntrl, only: LuTDM
      use Symmetry_Info, only: nSym=>nIrrep, MUL


      IMPLICIT NONE
#include "rassi.fh"
      CHARACTER(LEN=8) CHARTYPE!,LABEL
      Integer NSS
      Real*8 USOR(NSS,NSS),USOI(NSS,NSS)
      INTEGER ASS,BSS, iOpt
      Real*8 ROTMAT(3,3)
      Real*8 DENSOUT(6,nbst**2)

      Integer IOFF(8)
      Integer, Allocatable:: MAPST(:), MAPSP(:), MAPMS(:)
      Real*8, allocatable, target:: SDMXR(:), SDMXI(:),
     &                              SDMYR(:), SDMYI(:),
     &                              SDMZR(:), SDMZI(:)
      Real*8, allocatable, target:: SDMXR2(:), SDMXI2(:),
     &                              SDMYR2(:), SDMYI2(:),
     &                              SDMZR2(:), SDMZI2(:)
      Type A2_array
         Real*8, Pointer:: A2(:)
      End Type A2_Array
      Type (A2_array)  :: pZMR(3), pZMI(3)
      Type (A2_array)  :: pZMR2(3), pZMI2(3)
      Real*8, allocatable:: TMPR(:), TMPI(:), SCR(:), TDMZZ(:)

      REAL*8 CGY, CGX, CG0, TDM, S1, S2, SM1, SM2, FACT, CGM, CGP,
     &       URR, UIR, URL, UIL, DCLEBS
      Integer NBSTS, ITYPE, ISS, ISF, JOB, MPLET, MSPROJ, KSS, KSF,
     &        MPLETK, MSPROJK, LSS, LSF, MPLETL, MSPROJL, JOB1, JOB2,
     &        ISY12, IEMPTY, IDISK, IGO, ITD, NB1_I, NB1_F, ISY1,
     &        NB2_I, NB2_F, NB1, ISY2, ISY12_MA, NB2, J, I, IJ, IDIR

      nbsts=nbst**2

c VV: dummy initialization
      CGY=-1
      CGX=-1
      CG0=-1
c Initialize
      DENSOUT(:,:) = 0.0D00
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
      CALL mma_allocate(SDMXR,nbsts,Label='SDMXR')
      CALL mma_allocate(SDMXI,nbsts,Label='SDMXI')
      CALL mma_allocate(SDMYR,nbsts,Label='SDMYR')
      CALL mma_allocate(SDMYI,nbsts,Label='SDMYI')
      CALL mma_allocate(SDMZR,nbsts,Label='SDMZR')
      CALL mma_allocate(SDMZI,nbsts,Label='SDMZI')
      SDMXR(:)=0.0D0
      SDMXI(:)=0.0D0
      SDMYR(:)=0.0D0
      SDMYI(:)=0.0D0
      SDMZR(:)=0.0D0
      SDMZI(:)=0.0D0
      pZMR(1)%A2=>SDMXR(:)
      pZMR(2)%A2=>SDMYR(:)
      pZMR(3)%A2=>SDMZR(:)
      pZMI(1)%A2=>SDMXI(:)
      pZMI(2)%A2=>SDMYI(:)
      pZMI(3)%A2=>SDMZI(:)

      CALL mma_allocate(SDMXR2,nbsts,Label='SDMXR2')
      CALL mma_allocate(SDMXI2,nbsts,Label='SDMXI2')
      CALL mma_allocate(SDMYR2,nbsts,Label='SDMYR2')
      CALL mma_allocate(SDMYI2,nbsts,Label='SDMYI2')
      CALL mma_allocate(SDMZR2,nbsts,Label='SDMZR2')
      CALL mma_allocate(SDMZI2,nbsts,Label='SDMZI2')
      SDMXR2(:)=0.0D0
      SDMXI2(:)=0.0D0
      SDMYR2(:)=0.0D0
      SDMYI2(:)=0.0D0
      SDMZR2(:)=0.0D0
      SDMZI2(:)=0.0D0
      pZMR2(1)%A2=>SDMXR2(:)
      pZMR2(2)%A2=>SDMYR2(:)
      pZMR2(3)%A2=>SDMZR2(:)
      pZMI2(1)%A2=>SDMXI2(:)
      pZMI2(2)%A2=>SDMYI2(:)
      pZMI2(3)%A2=>SDMZI2(:)

      CALL mma_allocate(TMPR,nbsts,Label='TMPR')
      CALL mma_allocate(TMPI,nbsts,Label='TMPI')
      TMPR(:)=0.0D0
      TMPI(:)=0.0D0

      CALL mma_allocate(SCR,nbsts,Label='SCR')
c zeroed inside the loop

      CALL mma_allocate(TDMZZ,NTDMZZ,Label='TDMZZ')
      TDMZZ(:)=0.0D0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MAIN LOOP OVER KSF/LSF
C WRITTEN AS IN PRPROP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CORRESPONDING SPIN-FREE STATES OF THE
C REQUESTED SPIN STATES
c      ASF=MAPST(ASS)
c      BSF=MAPST(BSS)

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
        TDMZZ(:)=0.0D0


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
cccccccccccccc
c nbsq and NTDMZZ are the same, but they are the
c sizes of the symmetry-adapted matrices
c what we need is a NBST**2-sized matrix without symmetry
c Thus we expand the final TDM in C1 symmetry
c Leave the zero matrix elements as they are
ccccccccccccccc

cc Expand into C1
        ITD=0
        NB1_i=0
        NB1_f=0
        Do ISY1=1,NSYM
          NB2_i=0
          NB2_f=0
          NB1=NBASF(ISY1)
          NB1_f=NB1_i+NB1
          Do ISY2=1,NSYM
            ISY12_ma=MUL(ISY1,ISY2)
            NB2=NBASF(ISY2)
            NB2_f=NB2_i+NB2
            If(ISY12_ma.EQ.ISY12) then
              Do J=NB2_i+1,NB2_f
                Do I=NB1_i+1,NB1_f
                  ITD=ITD+1
                  TDM=TDMZZ(ITD)
                  IJ=I+NBST*(J-1)
                  SCR(IJ)=SCR(IJ)+TDM
                Enddo
              Enddo
            Endif
            NB2_i=NB2_i+NB2
          Enddo
          NB1_i=NB1_i+NB1
        Enddo

c ie, see how AMFI is processed in soeig.f
        SDMXR(:)=0.0D0
        SDMXI(:)=0.0D0
        SDMYR(:)=0.0D0
        SDMYI(:)=0.0D0
        SDMZR(:)=0.0D0
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
          CALL DAXPY_(nbsts,1.0d0,SCR,1,SDMZR,1)
        ELSE IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
          If (iOpt.eq.1) Then
          CALL DAXPY_(nbsts,CGX*ROTMAT(1,1),SCR,1,SDMXR,1)
          CALL DAXPY_(nbsts,CGY*ROTMAT(2,1),SCR,1,SDMXI,1)
          CALL DAXPY_(nbsts,CG0*ROTMAT(3,1),SCR,1,SDMXR,1)

          CALL DAXPY_(nbsts,CGX*ROTMAT(1,2),SCR,1,SDMYR,1)
          CALL DAXPY_(nbsts,CGY*ROTMAT(2,2),SCR,1,SDMYI,1)
          CALL DAXPY_(nbsts,CG0*ROTMAT(3,2),SCR,1,SDMYR,1)

          CALL DAXPY_(nbsts,CGX*ROTMAT(1,3),SCR,1,SDMZR,1)
          CALL DAXPY_(nbsts,CGY*ROTMAT(2,3),SCR,1,SDMZI,1)
          CALL DAXPY_(nbsts,CG0*ROTMAT(3,3),SCR,1,SDMZR,1)
          Else
          CALL DAXPY_(nbsts,CGX,SCR,1,SDMXR,1)
          CALL DAXPY_(nbsts,CGY,SCR,1,SDMYI,1)
          CALL DAXPY_(nbsts,CG0,SCR,1,SDMZR,1)
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
          TMPR(:)=0.0D0
          TMPI(:)=0.0D0

C right side
          CALL DAXPY_(nbsts,       URR,pZMR(IDIR)%A2,1,TMPR,1)
          CALL DAXPY_(nbsts,-1.0d0*UIR,pZMI(IDIR)%A2,1,TMPR,1)
          CALL DAXPY_(nbsts,       UIR,pZMR(IDIR)%A2,1,TMPI,1)
          CALL DAXPY_(nbsts,       URR,pZMI(IDIR)%A2,1,TMPI,1)

C left side
         CALL DAXPY_(nbsts,       URL,TMPR,1,pZMR2(IDIR)%A2,1)
         CALL DAXPY_(nbsts,       UIL,TMPI,1,pZMR2(IDIR)%A2,1)
         CALL DAXPY_(nbsts,       URL,TMPI,1,pZMI2(IDIR)%A2,1)
         CALL DAXPY_(nbsts,-1.0d0*UIL,TMPR,1,pZMI2(IDIR)%A2,1)
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
        DO I=1,nbsts
          DENSOUT(1,I)=SDMXR2(I)
          DENSOUT(2,I)=SDMYR2(I)
          DENSOUT(3,I)=SDMZR2(I)
          DENSOUT(4,I)=SDMXI2(I)
          DENSOUT(5,I)=SDMYI2(I)
          DENSOUT(6,I)=SDMZI2(I)
        END DO
      ELSE
        DO I=1,nbsts
          DENSOUT(1,I)=SDMZR2(I)
          DENSOUT(2,I)=SDMZR2(I)
          DENSOUT(3,I)=SDMZR2(I)
          DENSOUT(4,I)=SDMZI2(I)
          DENSOUT(5,I)=SDMZI2(I)
          DENSOUT(6,I)=SDMZI2(I)
        END DO
      END IF

c Free memory
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(TDMZZ)

      nullify(pZMI(3)%A2,pZMI(2)%A2,pZMI(1)%A2,pZMR(3)%A2,pZMR(2)%A2,
     &        pZMR(1)%A2)
      Call mma_deallocate(SDMZI)
      Call mma_deallocate(SDMZR)
      Call mma_deallocate(SDMYI)
      Call mma_deallocate(SDMYR)
      Call mma_deallocate(SDMXI)
      Call mma_deallocate(SDMXR)

      Call mma_deallocate(TMPI)
      Call mma_deallocate(TMPR)

      nullify(pZMI2(3)%A2,pZMI2(2)%A2,pZMI2(1)%A2,pZMR2(3)%A2,
     &        pZMR2(2)%A2,pZMR2(1)%A2)
      Call mma_deallocate(SDMZI2)
      Call mma_deallocate(SDMZR2)
      Call mma_deallocate(SDMYI2)
      Call mma_deallocate(SDMYR2)
      Call mma_deallocate(SDMXI2)
      Call mma_deallocate(SDMXR2)

      CALL mma_deallocate(MAPST)
      CALL mma_deallocate(MAPSP)
      CALL mma_deallocate(MAPMS)

      END SUBROUTINE MAKETDMAO
