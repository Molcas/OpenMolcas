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
      use constants, only: Zero
      use caspt2_global, only: FIMO, FAMO, FIFA, HONE, DREF
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM, NOMX, notri, NORB
      use definitions, only: iwp, wp, u6
      IMPLICIT None

      real(kind=wp), Allocatable:: BUF(:)
      integer(kind=iwp) IFTEST,NBUF,ISTLT,ISYM,NO

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
        WRITE(u6,*)'      TEST PRINTS FROM FOCK_RPT2.'
        WRITE(u6,*)'      ONE-ELECTRON HAMILTONIAN IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          IF ( NO.GT.0 ) THEN
            WRITE(u6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',HONE(ISTLT),NO)
            ISTLT=ISTLT+(NO*(NO+1))/2
          END IF
        END DO
      END IF

c Inactive and active Fock matrices:
      FIMO(:)=HONE(:)
      FAMO(:)=Zero
      CALL FMAT_CASPT2(FIMO,SIZE(FIMO),FAMO,SIZE(FAMO),DREF,SIZE(DREF),
     &                 NBUF,BUF)

* both FIMO and FAMO refer to the active space part only. FIMO comes
* from contractions over inactive orbitals, while FAMO from contractions
* over active orbitals and therefore are summed up together here
      FIFA(1:notri) = FIMO(1:notri)+FAMO(1:notri)

! these active orbital energies are not the ones used in
! MKFG3. Depending on whether the OUTO=canonical flag was set
! in &RASSCF, it will differ from the EPSA array in mkfg3.f

      IF ( IFTEST.NE.0 ) THEN
        WRITE(u6,*)'      INACTIVE FOCK MATRIX IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          IF ( NO.GT.0 ) THEN
            WRITE(u6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',FIMO(ISTLT),NO)
            ISTLT=ISTLT+(NO*(NO+1))/2
          END IF
        END DO

        WRITE(u6,*)'        ACTIVE FOCK MATRIX IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          IF ( NO.GT.0 ) THEN
            WRITE(u6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',FAMO(ISTLT),NO)
            ISTLT=ISTLT+(NO*(NO+1))/2
          END IF
        END DO

        WRITE(6,*)'      TOTAL FOCK MATRIX IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(u6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',FIFA(ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO

        ! these active orbital energies are not the ones used in
        ! MKFG3. Depending on whether the OUTO=canonical flag was set
        ! in &RASSCF, it will differ from the EPSA array in mkfg3.f
      END IF

      CALL mma_deallocate(BUF)

      END SUBROUTINE FOCK_RPT2
