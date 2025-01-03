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
      SUBROUTINE BINAT()
      use rassi_global_arrays, only : JBNUM, EIGVEC
      use rassi_aux, only : iDisk_TDM
      use OneDat, only: sNoNuc, sNoOri
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NBINA, NSTATE, IRREP, IBINA
      use cntrl, only: LuTDM
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      use rassi_data, only: NBSQ,NBMX,NBASF,NBST,NBTRI,NTDMZZ
      IMPLICIT NONE

      INTEGER IOFF_SEV, IOFF_VEC, IOFF_TDM, IOFF_ISV
      INTEGER IOPT, ICMP, ISYLAB, LS, LV, LV1
      INTEGER LE, I, LS1, ISEL, LS2, J, K, L, LE1
      INTEGER IJPAIR, KEIG_BRA, KEIG_KET, LSYM_BRA
      INTEGER LSYM_KET, LSYM12, IDISK, IV, IE, ITD, IRC, ISYM
      INTEGER ISYM1, ISYM2, NB, NB1, NB2, LV2, LE2, ITD1, ITD2
      INTEGER ISV, LB, LK, NBMIN, iEmpty, iGo
      INTEGER LUNIT, IDUMMY
      REAL*8  SSEL, SWAP, S_EV, X, DUMMY, SUMSNG
      CHARACTER(LEN=16) KNUM
      CHARACTER(LEN=8)  BNUM, LABEL
      CHARACTER(LEN=21) TXT
      CHARACTER(LEN=24) FNAME

C Tables of starting locations, created and used later
      DIMENSION IOFF_VEC(8),IOFF_SEV(8),IOFF_TDM(8),IOFF_ISV(8)
      DIMENSION IDUMMY(2),DUMMY(2)

      Integer, EXTERNAL :: ISFREEUNIT
      Real*8, EXTERNAL :: DDOT_

      Real*8, Allocatable:: ONBAS(:), SCR(:), SEV(:), SAO(:)
      Real*8, Allocatable:: UMAT(:), VTMAT(:), SVAL(:)
      Real*8, Allocatable:: SNGV1(:), SNGV2(:)
      Real*8, Allocatable:: TDMAT(:), TDMAO(:)
      Real*8, Allocatable:: BRABNO(:), KETBNO(:)

C Nr of basis functions, total
      NBSQ=0
      DO ISYM=1,NSYM
       NBSQ=NBSQ+NBASF(ISYM)**2
      END DO
C============================================================
C START BY CREATING A SET OF ORTHONORMAL VECTORS:
      CALL mma_allocate(ONBAS,NBSQ,Label='ONBAS')
C EIGENVALUES OF OVERLAP MATRIX:
      CALL mma_allocate(SEV,NBST,Label='SEV')
C READ ORBITAL OVERLAP MATRIX.
      CALL mma_allocate(SAO,NBTRI,Label='SAO')
      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL = 'MLTPL  0'
      CALL RDONE(IRC,IOPT,LABEL,ICMP,SAO,ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  BINAT ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF

C LOOP OVER SYMMETRY BLOCKS
C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=1
      LV=1
      LE=1
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        CALL DCOPY_(NB**2,[0.0D0],0,ONBAS(LV),1)
        CALL DCOPY_(NB,[1.0D0],0,ONBAS(LV),NB+1)
        CALL JACOB(SAO(LS),ONBAS(LV),NB,NB)
C SORT IN ORDER OF DECREASING EIGENVALUES.
        LS1=LS
        DO I=1,NB-1
         ISEL=I
         SSEL=SAO(LS1)
         LS2=LS1
         DO J=I+1,NB
          LS2=LS2+J
          IF(SAO(LS2).GT.SSEL) THEN
            ISEL=J
            SSEL=SAO(LS2)
          END IF
         END DO
         IF (ISEL.GT.I) THEN
          LS2=LS-1+(ISEL*(ISEL+1))/2
          SWAP=SAO(LS2)
          SAO(LS2)=SAO(LS1)
          SAO(LS1)=SWAP
          DO K=1,NB
           SWAP=ONBAS(LV-1+K+NB*(ISEL-1))
           ONBAS(LV-1+K+NB*(ISEL-1))=ONBAS(LV-1+K+NB*(I-1))
           ONBAS(LV-1+K+NB*(I-1))=SWAP
          END DO
         END IF
         LS1=LS1+I+1
        END DO
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO I=1,NB
          S_EV=SAO(LS1)
          SEV(LE1)=S_EV
          IF (S_EV.GT.1.0D-14) THEN
           X=1.0D00/SQRT(S_EV)
           CALL DSCAL_(NB,X,ONBAS(LV1),1)
          ELSE
           CALL DCOPY_(NB,[0.0D0],0,ONBAS(LV1),1)
          END IF
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
        END DO
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
      END DO
      CALL mma_deallocate(SAO)
C Starting at ONBAS there is now symmetry blocks of CMO arrays
C describing orthonormal vectors. In case the AO overlap matrix is
C (almost) singular, one or more vectors at the end of each symmetry
C block will be null vectors.
C============================================================

C Left and right singular vectors, used temporarily
C in calls to SVD routine. Also temporary, singular values.
      CALL mma_allocate(UMAT,NBMX**2,Label='UMAT')
      CALL mma_allocate(VTMAT,NBMX**2,Label='VTMAT')
      CALL mma_allocate(SVAL,NBMX,Label='SVAL')
C Final bra and ket singular value array:
      CALL mma_allocate(SNGV1,NBST,Label='SNGV1')
C An extra copy for singular values in a different order
C used until proper GV support for binatural orbitals.
      CALL mma_allocate(SNGV2,NBST,Label='SNGV2')
C The transition density matrix, symmetry-blocked
C Symmetry blocks may combine different symmetries, but the size
C is certainly less than or equal to NBSQ.
      CALL mma_allocate(TDMAT,NBSQ,Label='TDMAT')
C Same, read buffer
      CALL mma_allocate(TDMAO,NBSQ,Label='TDMAO')
C Temporary intermediate in matrix multiplies
C (Also used as temporary when transposing some TDMAO matrices)
      CALL mma_allocate(SCR,NBSQ,Label='SCR')
C The BRA and KET binatural orbitals:
      CALL mma_allocate(BRABNO,NBSQ,Label='BRABNO')
      CALL mma_allocate(KETBNO,NBSQ,Label='KETBNO')

C A long loop over eigenstate pairs:
      DO IJPAIR=1,NBINA
C Requested state pairs for computation: (OBSOLETE)
       KEIG_BRA=IBINA(1,IJPAIR)
       KEIG_KET=IBINA(2,IJPAIR)
C Get symmetries, via jobiph number for the states:
       LSYM_BRA=IRREP(JBNUM(KEIG_BRA))
       LSYM_KET=IRREP(JBNUM(KEIG_KET))
C Combined symmetry:
       LSYM12=MUL(LSYM_BRA,LSYM_KET)
C For relating left and right symmetry blocks, offset tables are
C needed for the singular values and for the TDM.
       ITD=0
       ISV=0
       DO ISYM1=1,NSYM
        IOFF_TDM(ISYM1)=ITD
        IOFF_ISV(ISYM1)=ISV
        ISYM2=MUL(ISYM1,LSYM12)
        ITD=ITD+NBASF(ISYM1)*NBASF(ISYM2)
        ISV=ISV+NBASF(ISYM1)
       END DO
       TDMAT(:)=0.0D0
C DOUBLE LOOP OVER RASSCF WAVE FUNCTIONS
       DO I=1,NSTATE
        IF (IRREP(JBNUM(I)).NE.LSYM_BRA) cycle
        DO J=1,NSTATE
         IF (IRREP(JBNUM(J)).NE.LSYM_KET) cycle
C PICK UP TRANSITION DENSITY MATRIX FOR THIS PAIR OF RASSCF STATES:
C WEIGHT WITH WHICH THEY CONTRIBUTE IS EIGVEC(I,KEIG_BRA)*EIGVEC(J,KEIG_KET).
         X=EIGVEC(KEIG_BRA,i)*EIGVEC(KEIG_KET,j)
         IDISK=iDisk_TDM(J,I,1)
         IEMPTY=iDisk_TDM(J,I,2)
         iOpt=2
         iGo=1
         If (IAND(iEMPTY,1).ne.0) Then
         IF (I.GT.J) THEN
            CALL dens2file(TDMAO,TDMAO,TDMAO,
     &                     nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
         ELSE
C Pick up conjugate TDM array, and transpose it into TDMAO.
            CALL dens2file(SCR,SCR,SCR,
     &                     nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
C Loop over the receiving side:
           DO ISYM1=1,NSYM
            ISYM2=MUL(ISYM1,LSYM12)
            NB1=NBASF(ISYM1)
            NB2=NBASF(ISYM2)
            DO K=1,NB1
             DO L=1,NB2
              TDMAO(IOFF_TDM(ISYM1)+K+NB1*(L-1))=
     &              SCR(IOFF_TDM(ISYM2)+L+NB2*(K-1))
             END DO
            END DO
           END DO
         END IF
         CALL DAXPY_(NBSQ,X,TDMAO,1,TDMAT,1)
         END IF
        END DO
       END DO
C TDMAT() now contains the transition density matrix in AO basis for
C the eigenstates.
C ------------------------------------------------------------------

C LOOP OVER SYMMETRY BLOCKS OF TDMAT.
C On the ket side, ISYM2 is not looping sequentially so we need
C tables of offsets:
       IV=0
       IE=0
       DO ISYM=1,NSYM
        IOFF_VEC(ISYM)=IV
        IOFF_SEV(ISYM)=IE
        IV=IV+NBASF(ISYM)**2
        IE=IE+NBASF(ISYM)
       END DO
       SNGV1(:)=0.0D0
       SNGV2(:)=0.0D0
       ITD=0
       DO ISYM1=1,NSYM
        ISYM2=MUL(ISYM1,LSYM12)
        NB1=NBASF(ISYM1)
        NB2=NBASF(ISYM2)
        LV1=1+IOFF_VEC(ISYM1)
        LV2=1+IOFF_VEC(ISYM2)
        LE1=1+IOFF_SEV(ISYM1)
        LE2=1+IOFF_SEV(ISYM2)
C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVECTORS OF THE OVERLAP MATRIX:
        CALL DGEMM_('N','N',NB1,NB2,NB2,1.0D0,
     &              TDMAT(1+ITD),NB1,ONBAS(LV2),NB2,
     &         0.0D0, SCR,NB1)
        CALL DGEMM_('T','N',NB1,NB2,NB1,1.0D0,
     &                ONBAS(LV1),NB1,SCR,NB1,
     &         0.0D0, TDMAT(1+ITD),NB1)
        ITD1=ITD
        DO I=1,NB1
         CALL DSCAL_(NB2,SEV(LE1-1+I),TDMAT(1+ITD1),NB1)
         ITD1=ITD1+1
        END DO
        ITD2=ITD
        DO I=1,NB2
         CALL DSCAL_(NB1,SEV(LE2-1+I),TDMAT(1+ITD2),1)
         ITD2=ITD2+NB1
        END DO

C SVD DECOMPOSITION OF THIS MATRIX BLOCK:
        CALL FULL_SVD(NB1,NB2,TDMAT(1+ITD),
     &                UMAT,VTMAT,SVAL)
* On return, UMAT has dimension (NB1,NB1),
* On return, VTMAT has dimension (NB2,NB2), transpose storage
        NBMIN=MIN(NB1,NB2)
C REEXPRESS THE SINGULAR VECTORS USING AO BASIS:
        LB=1+IOFF_VEC(ISYM1)
        LK=1+IOFF_VEC(ISYM2)
        CALL DGEMM_('N','N',NB1,NB1,NB1,1.0D0,
     &                ONBAS(LV1),NB1,UMAT,NB1,
     &          0.0D0,BRABNO(LB),NB1)
        CALL DGEMM_('N','T',NB2,NB2,NB2,1.0D0,
     &                ONBAS(LV2),NB2,VTMAT,NB2,
     &          0.0D0,KETBNO(LK),NB2)

C Move the singular values into their proper places:
        CALL DCOPY_(NBMIN,SVAL,1,SNGV1(1+IOFF_ISV(ISYM1)),1)
        CALL DCOPY_(NBMIN,SVAL,1,SNGV2(1+IOFF_ISV(ISYM2)),1)
        ITD=ITD+NB1*NB2
       END DO

C WRITE OUT THIS SET OF BI-NATURAL ORBITALS. THE FILES WILL BE NAMED
C BIORB.x_y, where x,y are KEIG_BRA and KEIG_KET.
C The BRA and KET orbitals will be written as alpha and beta, respectively,
C and the singular values will be written as "occupation numbers".

        WRITE(6,*)' Binatural singular values for the transition from'
        WRITE(6,*)' ket eigenstate KEIG_KET to bra eigenstate KEIG_BRA'
        WRITE(6,'(1x,I2,a,i2)') KEIG_BRA,' <-- ',KEIG_KET
        DO I=1,NSYM
          NB=NBASF(I)
          IF( NB.NE.0 ) THEN
            WRITE(6,'(A,I2)')' SYMMETRY SPECIES:',I
            LS=1+IOFF_SEV(I)
            WRITE(6,'(1X,10F8.5)')(SNGV1(LS-1+J),J=1,NB)
          ENDIF
        END DO

        WRITE(BNUM,'(I8)') KEIG_BRA
        BNUM=ADJUSTL(BNUM)
        WRITE(KNUM,'(I8)') KEIG_KET
        KNUM=ADJUSTL(KNUM)
        TXT=TRIM(BNUM)//' <-- '//TRIM(KNUM)
        KNUM=TRIM(BNUM)//'_'//TRIM(KNUM)

        FNAME='BIORB.'//KNUM
        WRITE(6,'(A,A)')' Orbitals are written onto file id = ',FNAME
        LUNIT=50
        LUNIT=ISFREEUNIT(LUNIT)
        CALL WRVEC_(FNAME,LUNIT,'CO',1,NSYM,NBASF,NBASF,
     &              BRABNO,KETBNO,SNGV1,SNGV2,
     &     DUMMY, DUMMY, IDUMMY,
     &     '* Binatural orbitals from transition '//TRIM(TXT), 0 )
        CLOSE(LUNIT)
        SUMSNG=DDOT_(SUM(NBASF),SNGV1,1,SNGV2,1)
        CALL ADD_INFO("BINAT",[SUMSNG],1,5)

C End of very long loop over eigenstate pairs.
      END DO

      WRITE(6,*) repeat('*',80)
      CALL mma_deallocate(ONBAS)
      CALL mma_deallocate(SEV)
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(UMAT)
      CALL mma_deallocate(VTMAT)
      CALL mma_deallocate(SVAL)
      CALL mma_deallocate(SNGV1)
      CALL mma_deallocate(SNGV2)
      CALL mma_deallocate(TDMAT)
      CALL mma_deallocate(TDMAO)
      CALL mma_deallocate(BRABNO)
      CALL mma_deallocate(KETBNO)

      END SUBROUTINE BINAT
