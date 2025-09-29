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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE FOCK_RPT2()
      use caspt2_global, only: FIMO, FAMO, FIFA, HONE, DREF
      use stdalloc, only: mma_allocate, mma_deallocate
      use ChoCASPT2
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"
#include "pt2_guga.fh"

      Real*8, Allocatable:: BUF(:)

c Purpose: Compute the standard Fock matrix which defines
c the PT2 orbitals and the standard H0 hamiltonian.
c Available input data are: Effective one-electron hamiltonian
c for non-frozen space, in MO basis, at HONE.
c Two-electron integrals in MO basis, second-order transformed,
c as the three integral sets on LUINTM.
c To be called from ORBCTL section, after second order two-el
c transformation, and TRAONE, are finished, or from H0CTL.

#ifdef _DEBUGPRINT_
      IFTEST=1
#else
      IFTEST=0
#endif

c notri=Size of an array with symmetry-blocked triangular
c submatrices, using non-frozen, non-deleted MO indices.
c NBUF=Max size of a LUINTM buffer.
      NBUF=MAX(NOMX**2,notri)
      CALL mma_allocate(BUF,NBUF,Label='BUF')

c One-electron Hamiltonian is in HONE

      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)'      TEST PRINTS FROM FOCK_RPT2.'
        WRITE(6,*)'      ONE-ELECTRON HAMILTONIAN IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          IF ( NO.GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',HONE(ISTLT),NO)
            ISTLT=ISTLT+(NO*(NO+1))/2
          END IF
        END DO
      END IF

c Inactive and active Fock matrices:
      FIMO(:)=HONE(:)
      FAMO(:)=0.0D0
      CALL FMAT_CASPT2(FIMO,SIZE(FIMO),FAMO,SIZE(FAMO),DREF,SIZE(DREF),
     &                 NBUF,BUF)

* both FIMO and FAMO refer to the active space part only. FIMO comes
* from contractions over inactive orbitals, while FAMO from contractions
* over active orbitals and therefore are summed up together here
      FIFA(1:notri) = FIMO(1:notri)+FAMO(1:notri)

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
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSI=IEPSI+1
          EPSI(IEPSI)=E
        END DO
        DO I=NI+1,NI+NA
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSA=IEPSA+1
          EPSA(IEPSA)=E
        END DO
        DO I=NI+NA+1,NO
          E=FIFA(ISTLT+(I*(I+1))/2)
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
      ! EASUM=0.0D0
      ! DO ISYM=1,NSYM
      !   NA=NASH(ISYM)
      !   DO I=1,NA
      !     ITOT=NAES(ISYM)+I
      !     ID=(ITOT*(ITOT+1))/2
      !     EASUM=EASUM+EPSA(ITOT)*DREF(ID)
      !   END DO
      ! END DO

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
        ISTLT=1
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',FIFA(ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO

        WRITE(6,*)
        WRITE(6,*)' FOCK_RPT2: ORBITAL ENERGIES, EPS:'
        WRITE(6,'(1X,5F12.6)')(EPS(I),I=1,NORBT)
        WRITE(6,*)'      INACTIVE ORBITAL ENERGIES, EPSI:'
        WRITE(6,'(1X,5F12.6)')(EPSI(I),I=1,NISHT)
        ! these active orbital energies are not the ones used in
        ! MKFG3. Depending on whether the OUTO=canonical flag was set
        ! in &RASSCF, it will differ from the EPSA array in mkfg3.f
        WRITE(6,*)'        ACTIVE ORBITAL ENERGIES, EPSA:'
        WRITE(6,'(1X,5F12.6)')(EPSA(I),I=1,NASHT)
        WRITE(6,*)'      EXTERNAL ORBITAL ENERGIES, EPSE:'
        WRITE(6,'(1X,5F12.6)')(EPSE(I),I=1,NSSHT)
      END IF

      CALL mma_deallocate(BUF)

      END SUBROUTINE FOCK_RPT2
