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
      SUBROUTINE NATSPIN_RASSI(DMAT,TDMZZ,VNAT,OCC,EIGVEC)
      use OneDat, only: sNoNuc, sNoOri
      use rassi_aux, only : iDisk_TDM
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: nState, NrNATO
      use cntrl, only: LuTDM
      use Symmetry_Info, only: nSym=>nIrrep

      IMPLICIT None
#include "rassi.fh"
      Real*8 DMAT(NBSQ),TDMZZ(NTDMZZ),VNAT(NBSQ),OCC(NBST)
      REAL*8 EIGVEC(NSTATE,NSTATE)

      CHARACTER(LEN=14) FNAME
      CHARACTER(LEN=8) KNUM, LABEL

      INTEGER, EXTERNAL :: ISFREEUNIT
      REAL*8, EXTERNAL :: DDOT_

      Real*8 Dummy(1)
      Integer iDummy(7,8)
      Real*8, allocatable:: SZZ(:), VEC(:), VEC2(:), SCR(:), EIG(:)

      Integer NSZZ, NVEC, NVEC2, NSCR, NEIG, IRC, IOPT, ICMP, ISYLAB,
     &        LS, LV, LE, ISYM, NB, LS1, LV1, LE1, I, KEIG, J, IEMPTY,
     &        IDISK, IGO, ID, INV, IOCC, ID1, ID2, ISCR, IJ, JI, I1,
     &        I2, ISTOCC, LUXXVEC, II
      REAL*8 X, SumOcc

C ALLOCATE WORKSPACE AREAS.
      NSZZ=NBTRI
      NVEC=NBSQ
      NVEC2=NBMX**2
      NSCR=NBMX**2
      NEIG=NBST
      CALL mma_allocate(SZZ,NSZZ,Label='SZZ')
      CALL mma_allocate(VEC,NVEC,Label='VEC')
      CALL mma_allocate(VEC2,NVEC2,Label='VEC2')
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      CALL mma_allocate(EIG,NEIG,Label='EIG')
C READ ORBITAL OVERLAP MATRIX.
      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      CALL RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  NATSPIN ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF
C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=1
      LV=1
      LE=1
      VEC(:)=0.0D0
      DO 100 ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO 20 I=1,NB**2,(NB+1)
          VEC(LV-1+I)=1.0D00
20      CONTINUE
        CALL JACOB(SZZ(LS),VEC(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO 30 I=1,NB
          EIG(LE1)=SZZ(LS1)
          X=1.0D00/SQRT(MAX(SZZ(LS1),1.0D-14))
          CALL DSCAL_(NB,X,VEC(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
30      CONTINUE
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
100   CONTINUE
      Call mma_deallocate(SZZ)

C VERY LONG LOOP OVER EIGENSTATES KEIG.
      DO KEIG=1,NRNATO

        CALL DCOPY_(NBSQ,[0.0D0],0,DMAT,1)
C DOUBLE LOOP OVER RASSCF WAVE FUNCTIONS, TRIANGULAR.
        DO I=1,NSTATE
          DO J=1,I
C WEIGHT WITH WHICH THIS TDM CONTRIBUTES IS EIGVEC(I,KEIG)*EIGVEC(J,KEIG).
C HOWEVER, WE ARE LOOPING TRIANGULARLY AND WILL RESTORE SYMMETRY BY
C ADDING TRANSPOSE AFTER DMAT HAS BEEN FINISHED, SO I=J IS SPECIAL CASE:
            X=EIGVEC(I,KEIG)*EIGVEC(J,KEIG)
            IF(ABS(X).GT.1.0D-12) THEN
              iDisk=iDisk_TDM(J,I,1)
              iEmpty=iDisk_TDM(J,I,2)
              If (IAND(iEmpty,2).ne.0) Then
                 iDisk=iDisk_TDM(J,I,1)
                 iOpt=2
                 iGo=2
C PICK UP TRANSITION SPIN DENSITY MATRIX FOR THIS PAIR OF RASSCF STATES:
                 CALL dens2file(TDMZZ,TDMZZ,TDMZZ,nTDMZZ,
     &                          LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
                 IF(I.EQ.J) X=0.5D00*X
                 CALL DAXPY_(NTDMZZ,X,TDMZZ,1,DMAT,1)
              End If
            END IF
          END DO
        END DO
C LOOP OVER SYMMETRY BLOCKS OF DMAT.
        ID=1
        INV=1
        IOCC=0
        LV=1
        LE=1
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVECTORS OF THE OVERLAP MATRIX:
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 DMAT(ID),NB,VEC(LV),NB,
     &           0.0D0,SCR,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,SCR,NB,
     &           0.0D0,DMAT(ID),NB)
          ID1=ID
          ID2=ID
          DO I=1,NB
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO
C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          ISCR=1
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
              SCR(ISCR)=DMAT(ID-1+IJ)+DMAT(ID-1+JI)
              ISCR=ISCR+1
            END DO
          END DO
C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          VEC2(:)=0.0D0
          CALL DCOPY_(NB,[1.0D0],0,VEC2,NB+1)
          CALL JACOB(SCR,VEC2,NB,NB)
          CALL JACORD(SCR,VEC2,NB,NB)
C JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=0
          DO I=1,NB
            II=II+I
            OCC(IOCC+NB+1-I)=SCR(II)
          END DO
          IOCC=IOCC+NB
C REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,VEC2,NB,
     &           0.0D0,SCR,NB)
          I1=1
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,SCR(I1),1,VNAT(I2),1)
            I1=I1+NB
          END DO
          ID=ID+NB**2
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
        END DO

C WRITE OUT THIS SET OF NATURAL SPIN ORBITALS. THE FILES WILL BE NAMED
C SSORB.1, SSORB.2, ...

       WRITE(KNUM,'(I8)') KEIG
       KNUM=ADJUSTL(KNUM)
       FNAME='SSORB.'//KNUM
       WRITE(6,'(A,I2)')' NATURAL SPIN ORBITALS FOR EIGENSTATE NR ',KEIG
       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ID = ',FNAME
       WRITE(6,'(A)')   ' OCCUPATION NUMBERS:'
        ISTOCC=0
        DO I=1,NSYM
          NB=NBASF(I)
          IF( NB.NE.0 ) THEN
            WRITE(6,'(A,I2)')' SYMMETRY SPECIES:',I
            WRITE(6,'(1X,10F8.5)')(OCC(ISTOCC+J),J=1,NB)
          ENDIF
          ISTOCC=ISTOCC+NB
        END DO
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)
        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &             VNAT, OCC, Dummy, iDummy,
     &     '* NATURAL SPIN ORBITALS FROM RASSI EIGENSTATE NR '//
     &     TRIM(KNUM))
        SUMOCC=DDOT_(SUM(NBASF),OCC,1,OCC,1)
        CALL ADD_INFO("NATSPIN",[SUMOCC],1,5)

C End of very long loop over eigenstates KEIG.
      END DO

      WRITE(6,*) repeat('*',80)
      CALL mma_deallocate(VEC)
      CALL mma_deallocate(VEC2)
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(EIG)
      END SUBROUTINE NATSPIN_RASSI
