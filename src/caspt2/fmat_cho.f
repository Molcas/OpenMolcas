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
      SUBROUTINE FMAT_CHO(CMO,NCMO,FFAO,FIAO,FAAO,HONE,NHONE,FIMO,NFIMO,
     &                                                       FAMO,NFAMO)
      use constants, only: Zero, One
      use caspt2_global, only: FIFA
      use caspt2_global, only: LUONEM
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NBTRI, IEOF1M, notri, NSYM, NBAS,
     &                         NORB, NFRO, IAD1M
      use definitions, only: iwp, wp
#ifdef _DEBUGPRINT_
      use definitions, only: u6
#endif
      IMPLICIT None
      integer(kind=iwp), intent(in):: NCMO, NHONE, NFIMO, NFAMO
      real(kind=wp), intent(in):: CMO(NCMO)
      real(kind=wp), intent(in):: FFAO(NBTRI),FIAO(NBTRI),FAAO(NBTRI)
      real(kind=wp), intent(out):: HONE(NHONE),FIMO(NFIMO),FAMO(NFAMO)

      real(kind=wp), allocatable:: SCR1(:), SCR2(:), SCR3(:)
      integer(kind=iwp) I, IDISK, IFAO, IJ, IOFMO, ISYM, J, LSC, LSCI,
     &                  NB, NBBMX, NBBT, NBOMX, NF, NO, NO_X, NOOMX
#ifdef _DEBUGPRINT_
      integer(kind=iwp) ISTLT
#endif

C THIS ROUTINE IS USED IF THE TWO-ELECTRON INTEGRALS ARE
C REPRESENTED BY CHOLESKY VECTORS:
C TRANSFORM FOCK MATRICES COMPUTED BY TRACHO
C TO MO BASIS FOR USE IN CASPT2.



      NBBT=0
      NBBMX=0
      NBOMX=0
      NOOMX=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       NO=NORB(ISYM)
       NBBT=NBBT+(NB*(NB+1))/2
       NBBMX=MAX(NBBMX,NB*NB)
       NBOMX=MAX(NBOMX,NB*NO)
       NOOMX=MAX(NOOMX,NO*NO)
      END DO

      CALL mma_allocate(SCR1,NBBMX,LABEL='SCR1')
      CALL mma_allocate(SCR2,NBOMX,LABEL='SCR2')
      CALL mma_allocate(SCR3,NOOMX,LABEL='SCR3')

      IFAO=1
      IOFMO=0
      LSC=1
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       If (NB.eq.0) CYCLE
       NO=NORB(ISYM)
       NO_X = Max(1,NO)
       NF=NFRO(ISYM)
       LSCI=LSC+NF*NB
* The frozen Fock matrix:
       CALL SQUARE(FFAO(IFAO),SCR1,NB,1,NB)
       CALL DGEMM_('N','N',NB,NO,NB, One,SCR1,NB,
     &            CMO(LSCI),NB,Zero,SCR2,NB)
       CALL DGEMM_('T','N',NO,NO,NB, One,CMO(LSCI),NB,
     &            SCR2,NB,Zero,SCR3,NO_X)
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         HONE(IOFMO+IJ)=SCR3(I+NO*(J-1))
        END DO
       END DO
* The inactive Fock matrix:
       CALL SQUARE(FIAO(IFAO),SCR1,NB,1,NB)
       CALL DGEMM_('N','N',NB,NO,NB, One,SCR1,NB,
     &            CMO(LSCI),NB,Zero,SCR2,NB)
       CALL DGEMM_('T','N',NO,NO,NB, One,CMO(LSCI),NB,
     &            SCR2,NB,Zero,SCR3,NO_X)
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         FIMO(IOFMO+IJ)=SCR3(I+NO*(J-1))
        END DO
       END DO
* The active Fock matrix:
       CALL SQUARE(FAAO(IFAO),SCR1,NB,1,NB)
       CALL DGEMM_('N','N',NB,NO,NB, One,SCR1,NB,
     &            CMO(LSCI),NB,Zero,SCR2,NB)
       CALL DGEMM_('T','N',NO,NO,NB, One,CMO(LSCI),NB,
     &            SCR2,NB,Zero,SCR3,NO_X)
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         FAMO(IOFMO+IJ)=SCR3(I+NO*(J-1))
        END DO
       END DO
       IFAO=IFAO+(NB*(NB+1))/2
       IOFMO=IOFMO+(NO*(NO+1))/2
       LSC=LSC+NB**2
      END DO

      CALL mma_deallocate(SCR1)
      CALL mma_deallocate(SCR2)
      CALL mma_deallocate(SCR3)

c Transformed frozen Fock matrix = Effective one-electron
* Hamiltonian HONE at IAD1M(3)
      IDISK=IEOF1M
      IAD1M(3)=IDISK
      CALL DDAFILE(LUONEM,1,HONE,notri,IDISK)
      IEOF1M=IDISK

      CALL DAXPY_(notri,One,HONE,1,FIMO,1)
      CALL DCOPY_(NOTRI,FIMO,1,FIFA,1)
      CALL DAXPY_(notri,One,FAMO,1,FIFA,1)

#ifdef _DEBUGPRINT_
        WRITE(6,*)'      INACTIVE FOCK MATRIX IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          IF ( NO.GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',FIMO(ISTLT),NO)
            ISTLT=ISTLT+(NO*(NO+1))/2
          END IF
        END DO

        WRITE(6,*)'        ACTIVE FOCK MATRIX IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          IF ( NO.GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',FAMO(ISTLT),NO)
            ISTLT=ISTLT+(NO*(NO+1))/2
          END IF
        END DO

        WRITE(6,*)'      TOTAL FOCK MATRIX IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',FIFA(ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO

#endif

      END SUBROUTINE FMAT_CHO
