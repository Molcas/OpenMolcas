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
      use rassi_aux, only : iDisk_TDM
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION DMAT(NBSQ),TDMZZ(NTDMZZ),VNAT(NBSQ),OCC(NBST)
      CHARACTER*14 FNAME
      CHARACTER*8 KNUM
      REAL*8 EIGVEC(NSTATE,NSTATE)

      EXTERNAL ISFREEUNIT
      EXTERNAL DDOT_

      DIMENSION Dummy(1),iDummy(7,8)

C ALLOCATE WORKSPACE AREAS.
      NSZZ=NBTRI
      NVEC=NBSQ
      NVEC2=NBMX**2
      NSCR=NBMX**2
      NEIG=NBST
      CALL GETMEM('SZZ   ','ALLO','REAL',LSZZ,NSZZ)
      CALL GETMEM('VEC   ','ALLO','REAL',LVEC,NVEC)
      CALL GETMEM('VEC2  ','ALLO','REAL',LVEC2,NVEC2)
      CALL GETMEM('SCR   ','ALLO','REAL',LSCR,NSCR)
      CALL GETMEM('EIG   ','ALLO','REAL',LEIG,NEIG)
C READ ORBITAL OVERLAP MATRIX.
      IRC=-1
      IOPT=6
      ICMP=1
      ISYLAB=1
      CALL RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSZZ),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  NATSPIN ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF
C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=LSZZ
      LV=LVEC
      LE=LEIG
      CALL FZERO(WORK(LVEC),NVEC)
      DO 100 ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO 20 I=1,NB**2,(NB+1)
          WORK(LV-1+I)=1.0D00
20      CONTINUE
        CALL JACOB(WORK(LS),WORK(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO 30 I=1,NB
          EIG=WORK(LS1)
          WORK(LE1)=EIG
          X=1.0D00/SQRT(MAX(EIG,1.0D-14))
          CALL DSCAL_(NB,X,WORK(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
30      CONTINUE
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
100   CONTINUE
      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NSZZ)

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
        LV=LVEC
        LE=LEIG
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVECTORS OF THE OVERLAP MATRIX:
*          CALL MXMA(DMAT(ID),1,NB,WORK(LV),1,NB,WORK(LSCR),1,NB,
*     *              NB,NB,NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 DMAT(ID),NB,WORK(LV),NB,
     &           0.0D0,WORK(LSCR),NB)
*          CALL MXMA(WORK(LV),NB,1,WORK(LSCR),1,NB,DMAT(ID),1,NB,
*     *              NB,NB,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LSCR),NB,
     &           0.0D0,DMAT(ID),NB)
          ID1=ID
          ID2=ID
          DO I=1,NB
            EIG=WORK(LE-1+I)
            CALL DSCAL_(NB,EIG,DMAT(ID1),NB)
            CALL DSCAL_(NB,EIG,DMAT(ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO
C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          ISCR=LSCR
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
              WORK(ISCR)=DMAT(ID-1+IJ)+DMAT(ID-1+JI)
              ISCR=ISCR+1
            END DO
          END DO
C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          CALL DCOPY_(NVEC2,[0.0D0],0,WORK(LVEC2),1)
          CALL DCOPY_(NB,[1.0D0],0,WORK(LVEC2),NB+1)
          CALL JACOB(WORK(LSCR),WORK(LVEC2),NB,NB)
          CALL JACORD(WORK(LSCR),WORK(LVEC2),NB,NB)
C JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=LSCR-1
          DO I=1,NB
            II=II+I
            OCC(IOCC+NB+1-I)=WORK(II)
          END DO
          IOCC=IOCC+NB
C REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
*          CALL MXMA(WORK(LV),1,NB,WORK(LVEC2),1,NB,
*     *              WORK(LSCR),1,NB,NB,NB,NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LVEC2),NB,
     &           0.0D0,WORK(LSCR),NB)
          I1=LSCR
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,WORK(I1),1,VNAT(I2),1)
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
        IFOCC=1
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)
        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &     VNAT, OCC, Dummy, iDummy,
     &     '* NATURAL SPIN ORBITALS FROM RASSI EIGENSTATE NR '//
     &     TRIM(KNUM))
        SUMOCC=DDOT_(SUM(NBASF),OCC,1,OCC,1)
        CALL ADD_INFO("NATSPIN",[SUMOCC],1,5)

C End of very long loop over eigenstates KEIG.
      END DO

      WRITE(6,*)('*',I=1,80)
      CALL GETMEM('VEC   ','FREE','REAL',LVEC,NVEC)
      CALL GETMEM('VEC2  ','FREE','REAL',LVEC2,NVEC2)
      CALL GETMEM('SCR   ','FREE','REAL',LSCR,NSCR)
      CALL GETMEM('EIG   ','FREE','REAL',LEIG,NEIG)
      RETURN
      END
