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
      SUBROUTINE FMAT_CHO(CMO,FFAO,FIAO,FAAO,HONE,FIMO,FAMO)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "caspt2.fh"
#include "output.fh"
#include "SysDef.fh"
      DIMENSION  CMO(NCMO)
      DIMENSION FFAO(NBTRI),FIAO(NBTRI),FAAO(NBTRI)
      DIMENSION HONE(NHONE),FIMO(NFIMO),FAMO(NFAMO)

C THIS ROUTINE IS USED IF THE TWO-ELECTRON INTEGRALS ARE
C REPRESENTED BY CHOLESKY VECTORS:
C TRANSFORM FOCK MATRICES COMPUTED BY TRACHO
C TO MO BASIS FOR USE IN CASPT2.

      CALL QENTER('FMAT_CHO')

#ifdef _DEBUG_
      IFTEST=1
#else
      IFTEST=0
#endif

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

      CALL GETMEM('SCR1','Allocate','Real',LSCR1,NBBMX)
      CALL GETMEM('SCR2','Allocate','Real',LSCR2,NBOMX)
      CALL GETMEM('SCR3','Allocate','Real',LSCR3,NOOMX)

      IFAO=1
      IOFMO=0
      LSC=1
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       If (NB.eq.0) Go To 99
       NO=NORB(ISYM)
       NO_X = Max(1,NO)
       NF=NFRO(ISYM)
       LSCI=LSC+NF*NB
* The frozen Fock matrix:
       CALL SQUARE(FFAO(IFAO),WORK(LSCR1),NB,1,NB)
       CALL DGEMM_('N','N',NB,NO,NB, 1.0D0,WORK(LSCR1),NB,
     &            CMO(LSCI),NB,0.0D0,WORK(LSCR2),NB)
       CALL DGEMM_('T','N',NO,NO,NB, 1.0D0,CMO(LSCI),NB,
     &            WORK(LSCR2),NB,0.0D0,WORK(LSCR3),NO_X)
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         HONE(IOFMO+IJ)=WORK(LSCR3+I-1+NO*(J-1))
        END DO
       END DO
* The inactive Fock matrix:
       CALL SQUARE(FIAO(IFAO),WORK(LSCR1),NB,1,NB)
       CALL DGEMM_('N','N',NB,NO,NB, 1.0D0,WORK(LSCR1),NB,
     &            CMO(LSCI),NB,0.0D0,WORK(LSCR2),NB)
       CALL DGEMM_('T','N',NO,NO,NB, 1.0D0,CMO(LSCI),NB,
     &            WORK(LSCR2),NB,0.0D0,WORK(LSCR3),NO_X)
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         FIMO(IOFMO+IJ)=WORK(LSCR3+I-1+NO*(J-1))
        END DO
       END DO
* The active Fock matrix:
       CALL SQUARE(FAAO(IFAO),WORK(LSCR1),NB,1,NB)
       CALL DGEMM_('N','N',NB,NO,NB, 1.0D0,WORK(LSCR1),NB,
     &            CMO(LSCI),NB,0.0D0,WORK(LSCR2),NB)
       CALL DGEMM_('T','N',NO,NO,NB, 1.0D0,CMO(LSCI),NB,
     &            WORK(LSCR2),NB,0.0D0,WORK(LSCR3),NO_X)
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         FAMO(IOFMO+IJ)=WORK(LSCR3+I-1+NO*(J-1))
        END DO
       END DO
       IFAO=IFAO+(NB*(NB+1))/2
       IOFMO=IOFMO+(NO*(NO+1))/2
       LSC=LSC+NB**2
 99    CONTINUE
      END DO

      CALL GETMEM('SCR1','Free','Real',LSCR1,NBBMX)
      CALL GETMEM('SCR2','Free','Real',LSCR2,NBOMX)
      CALL GETMEM('SCR3','Free','Real',LSCR3,NOOMX)

c Transformed frozen Fock matrix = Effective one-electron
* Hamiltonian HONE at IAD1M(3)
      IDISK=IEOF1M
      IAD1M(3)=IDISK
      CALL DDAFILE(LUONEM,1,HONE,notri,IDISK)
      IEOF1M=IDISK

      CALL DAXPY_(notri,1.0D00,HONE,1,FIMO,1)
      CALL DCOPY_(NOTRI,FIMO,1,WORK(LFIFA),1)
      CALL DAXPY_(notri,1.0D00,FAMO,1,WORK(LFIFA),1)

c   Orbital energies, EPS, EPSI,EPSA,EPSE:
      IEPS=0
      IEPSI=0
      IEPSA=0
      IEPSE=0
      ISTLT=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NO=NORB(ISYM)
        DO I=1,NI
          E=WORK(LFIFA+ISTLT-1+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSI=IEPSI+1
          EPSI(IEPSI)=E
        END DO
        DO I=NI+1,NI+NA
          E=WORK(LFIFA+ISTLT-1+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSA=IEPSA+1
          EPSA(IEPSA)=E
        END DO
        DO I=NI+NA+1,NO
          E=WORK(LFIFA+ISTLT-1+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSE=IEPSE+1
          EPSE(IEPSE)=E
        END DO
        ISTLT=ISTLT+(NO*(NO+1))/2
      END DO

C EASUM=CONTRACT EPSA WITH DIAGONAL OF ACTIVE DENS
C This is never used anywhere, and it is actually
C wrong in XMS, since the DREF used is not the average
C density.
      EASUM=0.0D0
      DO ISYM=1,NSYM
        NA=NASH(ISYM)
        DO I=1,NA
          ITOT=NAES(ISYM)+I
          ID=(ITOT*(ITOT+1))/2
          EASUM=EASUM+EPSA(ITOT)*WORK(LDREF-1+ID)
        END DO
      END DO

      IF ( IFTEST.NE.0 ) THEN
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
        ISTLT=0
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',WORK(LFIFA+ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO

        WRITE(6,*)
        WRITE(6,*)'      ORBITAL ENERGIES, EPS:'
        WRITE(6,'(1X,5F12.6)')(EPS(I),I=1,NORBT)
        WRITE(6,*)'      INACTIVE ORBITAL ENERGIES, EPSI:'
        WRITE(6,'(1X,5F12.6)')(EPSI(I),I=1,NISHT)
        WRITE(6,*)'        ACTIVE ORBITAL ENERGIES, EPSA:'
        WRITE(6,'(1X,5F12.6)')(EPSA(I),I=1,NASHT)
        WRITE(6,*)'      EXTERNAL ORBITAL ENERGIES, EPSE:'
        WRITE(6,'(1X,5F12.6)')(EPSE(I),I=1,NSSHT)
      END IF

      CALL QEXIT('FMAT_CHO')

      RETURN
      END
