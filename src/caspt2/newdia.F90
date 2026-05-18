!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
      SUBROUTINE NEWDIA()
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nSym, nInDep, nASup, nISup
      use constants, only: Zero
      use definitions, only: iwp, wp
      IMPLICIT NONE

      INTEGER(kind=iwp) ICASE,ISYM,NIN,NAS,NIS,I
      INTEGER(kind=iwp) JD
      REAL(kind=wp), ALLOCATABLE:: BD(:), ID(:), C1(:), C2(:)

! Post-diagonalization modification of diagonal energy
! denominator terms for active and for non-active superindex.


      DO ICASE=1,13
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) CYCLE
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) CYCLE
! Remember: NIN values in BDIAG, but must read NAS for correct
! positioning.
          CALL mma_allocate(BD,NAS,LABEL='BD')
          CALL mma_allocate(ID,NIS,LABEL='ID')
          CALL mma_allocate(C1,NAS,LABEL='C1')
          CALL mma_allocate(C2,NIS,LABEL='C2')
          JD=IDBMAT(ISYM,ICASE)
! Active, and non-active, energy denominators:
          CALL DDAFILE(LUSBT,2,BD,NAS,JD)
          CALL DDAFILE(LUSBT,2,ID,NIS,JD)
! Active, and non-active, corrections:
! (Replace this strange example with something sensible)
          C1(:)=Zero
          C2(:)=Zero
! Modifications are added to the usual diagonal energies:
          DO I=1,NAS
            BD(I)=BD(I)+C1(I)
          END DO
          DO I=1,NIS
            ID(I)=ID(I)+C2(I)
          END DO
          JD=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,BD,NAS,JD)
          CALL DDAFILE(LUSBT,1,ID,NIS,JD)
! Added modifications are saved on LUSBT.
          CALL DDAFILE(LUSBT,1,C1,NAS,JD)
          CALL DDAFILE(LUSBT,1,C2,NIS,JD)

          CALL mma_deallocate(BD)
          CALL mma_deallocate(ID)
          CALL mma_deallocate(C1)
          CALL mma_deallocate(C2)
          JD=IDBMAT(ISYM,ICASE)
        END DO
      END DO
      END SUBROUTINE NEWDIA
