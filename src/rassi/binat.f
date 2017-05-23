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
      IMPLICIT NONE

#include "SysDef.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      INTEGER IOFF_SEV, IOFF_VEC, IOFF_TDM, IOFF_ISV
      INTEGER LONBAS, LSEV, LSAO, IOPT, ICMP, ISYLAB, LS, LV
      INTEGER LE, I, LS1, ISEL, LS2, J, K, L, LV1, LE1
      INTEGER LUMAT, LVTMAT, LSVAL, LTDMAT, LTDMAO, LSCR
      INTEGER LBRABNO, LKETBNO, LSNGV1, LSNGV2
      INTEGER IJPAIR, KEIG_BRA, KEIG_KET, LSYM_BRA
      INTEGER LSYM_KET, LSYM12, IDISK, IV, IE, ITD, IRC, ISYM
      INTEGER ISYM1, ISYM2, NB, NB1, NB2, LV2, LE2, ITD1, ITD2
      INTEGER ISV, LB, LK, NBMIN
      INTEGER LUNIT, ISFREEUNIT, IDUMMY
      REAL*8  SSEL, SWAP, SEV, X, DUMMY
      CHARACTER*4  KNUM
      CHARACTER*9  TXT
      CHARACTER*10 FNAME

C Tables of starting locations, created and used later
      DIMENSION IOFF_VEC(8),IOFF_SEV(8),IOFF_TDM(8),IOFF_ISV(8)
      DIMENSION IDUMMY(2),DUMMY(2)

      EXTERNAL ISFREEUNIT

      Call qEnter('BINAT')
C Nr of basis functions, total
      NBSQ=0
      DO ISYM=1,NSYM
       NBSQ=NBSQ+NBASF(ISYM)**2
      END DO
C============================================================
C START BY CREATING A SET OF ORTHONORMAL VECTORS:
      CALL GETMEM('ONBAS','ALLO','REAL',LONBAS,NBSQ)
C EIGENVALUES OF OVERLAP MATRIX:
      CALL GETMEM('SEV','ALLO','REAL',LSEV,NBST)
C READ ORBITAL OVERLAP MATRIX.
      CALL GETMEM('SAO','ALLO','REAL',LSAO,NBTRI)
      IRC=-1
      IOPT=6
      ICMP=1
      ISYLAB=1
      CALL RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSAO),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  BINAT ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF

C LOOP OVER SYMMETRY BLOCKS
C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=LSAO
      LV=LONBAS
      LE=LSEV
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        CALL DCOPY_(NB**2,0.0D0,0,WORK(LV),1)
        CALL DCOPY_(NB,1.0D0,0,WORK(LV),NB+1)
        CALL JACOB(WORK(LS),WORK(LV),NB,NB)
C SORT IN ORDER OF DECREASING EIGENVALUES.
        LS1=LS
        DO I=1,NB-1
         ISEL=I
         SSEL=WORK(LS1)
         LS2=LS1
         DO J=I+1,NB
          LS2=LS2+J
          IF(WORK(LS2).GT.SSEL) THEN
            ISEL=J
            SSEL=WORK(LS2)
          END IF
         END DO
         IF (ISEL.GT.I) THEN
          LS2=LS-1+(ISEL*(ISEL+1))/2
          SWAP=WORK(LS2)
          WORK(LS2)=WORK(LS1)
          WORK(LS1)=SWAP
          DO K=1,NB
           SWAP=WORK(LV-1+K+NB*(ISEL-1))
           WORK(LV-1+K+NB*(ISEL-1))=WORK(LV-1+K+NB*(I-1))
           WORK(LV-1+K+NB*(I-1))=SWAP
          END DO
         END IF
         LS1=LS1+I+1
        END DO
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO I=1,NB
          SEV=WORK(LS1)
          WORK(LE1)=SEV
          IF (SEV.GT.1.0D-14) THEN
           X=1.0D00/SQRT(SEV)
           CALL DSCAL_(NB,X,WORK(LV1),1)
          ELSE
           CALL DCOPY_(NB,0.0D0,0,WORK(LV1),1)
          END IF
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
        END DO
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
      END DO
      CALL GETMEM('SAO','FREE','REAL',LSAO,NBTRI)
C Starting at Work(LONBAS) there is now symmetry blocks of CMO arrays
C describing orthonormal vectors. In case the AO overlap matrix is
C (almost) singular, one or more vectors at the end of each symmetry
C block will be null vectors.
C============================================================

C Left and right singular vectors, used temporarily
C in calls to SVD routine. Also temporary, singular values.
      CALL GETMEM('UMAT','ALLO','REAL',LUMAT,NBMX**2)
      CALL GETMEM('VTMAT','ALLO','REAL',LVTMAT,NBMX**2)
      CALL GETMEM('SVEC','ALLO','REAL',LSVAL,NBMX)
C Final bra and ket singular value array:
      CALL GETMEM('SNGV1','ALLO','REAL',LSNGV1,NBST)
C An extra copy for singular values in a different order
C used until proper GV support for binatural orbitals.
      CALL GETMEM('SNGV2','ALLO','REAL',LSNGV2,NBST)
C The transition density matrix, symmetry-blocked
C Symmetry blocks may combine different symmetries, but the size
C is certainly less than or equal to NBSQ.
      CALL GETMEM('TDMAT','ALLO','REAL',LTDMAT,NBSQ)
C Same, read buffer
      CALL GETMEM('TDMAO','ALLO','REAL',LTDMAO,NBSQ)
C Temporary intermediate in matrix multiplies
C (Also used as temporary when transposing some TDMAO matrices)
      CALL GETMEM('SCR','ALLO','REAL',LSCR,NBSQ)
C The BRA and KET binatural orbitals:
      CALL GETMEM('BRABNO','ALLO','REAL',LBRABNO,NBSQ)
      CALL GETMEM('KETBNO','ALLO','REAL',LKETBNO,NBSQ)

C A long loop over eigenstate pairs:
      DO IJPAIR=1,NBINA
C Requested state pairs for computation:
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
       CALL DCOPY_(NBSQ,0.0D0,0,WORK(LTDMAT),1)
C DOUBLE LOOP OVER RASSCF WAVE FUNCTIONS
       DO I=1,NSTATE
        IF (IRREP(JBNUM(I)).NE.LSYM_BRA) GOTO 92
        DO J=1,NSTATE
         IF (IRREP(JBNUM(J)).NE.LSYM_KET) GOTO 91
C PICK UP TRANSITION DENSITY MATRIX FOR THIS PAIR OF RASSCF STATES:
C WEIGHT WITH WHICH THEY CONTRIBUTE IS EIGVEC(I,KEIG_BRA)*EIGVEC(J,KEIG_KET).
         X=EIGVEC(I,KEIG_BRA)*EIGVEC(J,KEIG_KET)
         IF (I.GT.J) THEN
           IDISK=IDTDM(I,J)
           CALL DDAFILE(LUTDM,2,WORK(LTDMAO),NBSQ,IDISK)
         ELSE
C Pick up conjugate TDM array, and transpose it into TDMAO.
           IDISK=IDTDM(J,I)
           CALL DDAFILE(LUTDM,2,WORK(LSCR),NBSQ,IDISK)
C Loop over the receiving side:
           DO ISYM1=1,NSYM
            ISYM2=MUL(ISYM1,LSYM12)
            NB1=NBASF(ISYM1)
            NB2=NBASF(ISYM2)
            DO K=1,NB1
             DO L=1,NB2
              WORK(LTDMAO-1+IOFF_TDM(ISYM1)+K+NB1*(L-1))=
     &              WORK(LSCR-1+IOFF_TDM(ISYM2)+L+NB2*(K-1))
             END DO
            END DO
           END DO
         END IF
         CALL DAXPY_(NBSQ,X,WORK(LTDMAO),1,WORK(LTDMAT),1)
   91    CONTINUE
        END DO
   92   CONTINUE
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
       CALL DCOPY_(NBST,0.0D0,0,WORK(LSNGV1),1)
       CALL DCOPY_(NBST,0.0D0,0,WORK(LSNGV2),1)
       ITD=0
       DO ISYM1=1,NSYM
        ISYM2=MUL(ISYM1,LSYM12)
        NB1=NBASF(ISYM1)
        NB2=NBASF(ISYM2)
        LV1=LONBAS+IOFF_VEC(ISYM1)
        LV2=LONBAS+IOFF_VEC(ISYM2)
        LE1=LSEV+IOFF_SEV(ISYM1)
        LE2=LSEV+IOFF_SEV(ISYM2)
C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVECTORS OF THE OVERLAP MATRIX:
        CALL DGEMM_('N','N',NB1,NB2,NB2,1.0D0,
     &              WORK(LTDMAT+ITD),NB1,WORK(LV2),NB2,
     &         0.0D0, WORK(LSCR),NB1)
        CALL DGEMM_('T','N',NB1,NB2,NB1,1.0D0,
     &                WORK(LV1),NB1,WORK(LSCR),NB1,
     &         0.0D0, WORK(LTDMAT+ITD),NB1)
        ITD1=ITD
        DO I=1,NB1
         CALL DSCAL_(NB2,WORK(LE1-1+I),WORK(LTDMAT+ITD1),NB1)
         ITD1=ITD1+1
        END DO
        ITD2=ITD
        DO I=1,NB2
         CALL DSCAL_(NB1,WORK(LE2-1+I),WORK(LTDMAT+ITD2),1)
         ITD2=ITD2+NB1
        END DO

C SVD DECOMPOSITION OF THIS MATRIX BLOCK:
        CALL FULL_SVD(NB1,NB2,WORK(LTDMAT+ITD),
     &                WORK(LUMAT),WORK(LVTMAT),WORK(LSVAL))
* On return, WORK(LUMAT) has dimension (NB1,NB1),
* On return, WORK(LVTMAT) has dimension (NB2,NB2), transpose storage
        NBMIN=MIN(NB1,NB2)
C REEXPRESS THE SINGULAR VECTORS USING AO BASIS:
        LB=LBRABNO+IOFF_VEC(ISYM1)
        LK=LKETBNO+IOFF_VEC(ISYM2)
        CALL DGEMM_('N','N',NB1,NB1,NB1,1.0D0,
     &                WORK(LV1),NB1,WORK(LUMAT),NB1,
     &          0.0D0,WORK(LB),NB1)
        CALL DGEMM_('N','T',NB2,NB2,NB2,1.0D0,
     &                WORK(LV2),NB2,WORK(LVTMAT),NB2,
     &          0.0D0,WORK(LK),NB2)

C Move the singular values into their proper places:
        CALL DCOPY_(NBMIN,WORK(LSVAL),1,WORK(LSNGV1+IOFF_ISV(ISYM1)),1)
        CALL DCOPY_(NBMIN,WORK(LSVAL),1,WORK(LSNGV2+IOFF_ISV(ISYM2)),1)
        ITD=ITD+NB1*NB2
       END DO

C WRITE OUT THIS SET OF BI-NATURAL ORBITALS. THE FILES WILL BE NAMED
C BRAORB0x0y, where x,y are KEIG_BRA and KEIG_KET, and KETORB0x0y,
C similarly. Later, when there is GV support for binatural orbital
C pairs, there will be one single orbital file.
C Using separate conventional files, the singular values will be
C written as ''occupation numbers''.

        WRITE(6,*)' Binatural singular values for the transition from'
        WRITE(6,*)' ket eigenstate KEIG_KET to bra eigenstate KEIG_BRA'
        WRITE(6,'(1x,I2,a,i2)') KEIG_BRA,' <-- ',KEIG_KET
        DO I=1,NSYM
          NB=NBASF(I)
          IF( NB.NE.0 ) THEN
            WRITE(6,'(A,I2)')' SYMMETRY SPECIES:',I
            LS=LSNGV1+IOFF_SEV(I)
            WRITE(6,'(1X,10F8.5)')(WORK(LS-1+J),J=1,NB)
          ENDIF
        END DO

        WRITE(6,*)' Presently, GV lacks support for viewing binatural'
        WRITE(6,*)' pairs of orbitals from a single orbital files.'
        WRITE(6,*)' For the moment, two individual files are produced:'
        WRITE(KNUM,'(I2.2,I2.2)') KEIG_BRA,KEIG_KET
        TXT=KNUM(1:2)//' <-- '//KNUM(3:4)

        FNAME='BRAORB'//KNUM
        WRITE(6,*)' File ',FNAME,' with the left singular orbitals,'
        LUNIT=50
        LUNIT=ISFREEUNIT(LUNIT)
        CALL WRVEC(FNAME,LUNIT,'CO',NSYM,NBASF,NBASF,
     &     WORK(LBRABNO),WORK(LSNGV1), DUMMY, IDUMMY,
     &     '* Left singular orbitals from transition '//TXT )
        CLOSE(LUNIT)

        FNAME='KETORB'//KNUM
        WRITE(6,*)' and ',FNAME,' with the right singular orbitals.'
        LUNIT=50
        LUNIT=ISFREEUNIT(LUNIT)
        CALL WRVEC(FNAME,LUNIT,'CO',NSYM,NBASF,NBASF,
     &     WORK(LKETBNO),WORK(LSNGV2), DUMMY, IDUMMY,
     &     '* Right singular orbitals from transition '//TXT )
        CLOSE(LUNIT)
C End of very long loop over eigenstate pairs.
      END DO

      WRITE(6,*)('*',I=1,80)
      CALL GETMEM('ONBAS','FREE','REAL',LONBAS,NBSQ)
      CALL GETMEM('SEV   ','FREE','REAL',LSEV,NBST)
      CALL GETMEM('SCR','FREE','REAL',LSCR,NBSQ)
      CALL GETMEM('UMAT','FREE','REAL',LUMAT,NBMX**2)
      CALL GETMEM('VTMAT','FREE','REAL',LVTMAT,NBMX**2)
      CALL GETMEM('SVAL','FREE','REAL',LSVAL,NBMX)
      CALL GETMEM('SNGV1','FREE','REAL',LSNGV1,NBST)
      CALL GETMEM('SNGV2','FREE','REAL',LSNGV2,NBST)
      CALL GETMEM('TDMAT','FREE','REAL',LTDMAT,NBSQ)
      CALL GETMEM('TDMAO','FREE','REAL',LTDMAO,NBSQ)
      CALL GETMEM('BRABNO','FREE','REAL',LBRABNO,NBSQ)
      CALL GETMEM('KETBNO','FREE','REAL',LKETBNO,NBSQ)
      Call qExit('BINAT')
      RETURN
      END
