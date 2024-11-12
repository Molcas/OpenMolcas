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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
************************************************************************
* This file contains boiler-plate code to allow developers to experiment
* with modifications to the zero-order hamiltonian. To enable them,
* remove the #if 0/#endif blocks, and use 'HZero = Custom' in the input.
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE NEWB()
#if 0
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
#include "caspt2.fh"

      INTEGER ICASE,ISYM,NAS,NIS,NCOEF
      INTEGER IDS,NS,IDB,NB
      REAL*8, ALLOCATABLE:: S(:), B(:)

C Modify B matrices, if requested.


      DO ICASE=1,11
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NCOEF=NAS*NIS
          IF(NCOEF.EQ.0) CYCLE
          NS=(NAS*(NAS+1))/2
          CALL mma_allocate(S,NS,LABEL='S')
          IDS=IDSMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,S,NS,IDS)
          NB=NS
          CALL mma_allocate(B,NB,LABEL='B')
          IDB=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,B,NB,IDB)
C Modify B matrix, using S matrix and some other data.
          CALL mma_deallocate(B)
          CALL mma_deallocate(S)
        END DO
      END DO

#endif
      END SUBROUTINE NEWB

      SUBROUTINE NEWDIA()
#if 0
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
#include "caspt2.fh"

      INTEGER ICASE,ISYM,NIN,NAS,NIS,I
      INTEGER JD
      REAL*8, ALLOCATABLE:: BD(:), ID(:), C1(:), C2(:)

C Post-diagonalization modification of diagonal energy
C denominator terms for active and for non-active superindex.


      DO ICASE=1,13
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) CYCLE
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) CYCLE
C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          CALL mma_allocate(BD,NAS,LABEL='BD')
          CALL mma_allocate(ID,NIS,LABEL='ID')
          CALL mma_allocate(C1,NAS,LABEL='C1')
          CALL mma_allocate(C2,NIS,LABEL='C2')
          JD=IDBMAT(ISYM,ICASE)
C Active, and non-active, energy denominators:
          CALL DDAFILE(LUSBT,2,BD,NAS,JD)
          CALL DDAFILE(LUSBT,2,ID,NIS,JD)
C Active, and non-active, corrections:
C (Replace this strange example with something sensible)
          C1(:)=0.0D0
          C2(:)=0.0D0
C Modifications are added to the usual diagonal energies:
          DO I=1,NAS
            BD(I)=BD(I)+C1(I)
          END DO
          DO I=1,NIS
            ID(I)=ID(I)+C2(I)
          END DO
          ID=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,BD,NAS,JD)
          CALL DDAFILE(LUSBT,1,ID,NIS,JD)
C Added modifications are saved on LUSBT.
          CALL DDAFILE(LUSBT,1,C1,NAS,JD)
          CALL DDAFILE(LUSBT,1,C2,NIS,JD)

          CALL mma_deallocate(BD)
          CALL mma_deallocate(ID)
          CALL mma_deallocate(C1)
          CALL mma_deallocate(C2)
          ID=IDBMAT(ISYM,ICASE)
        END DO
      END DO
#endif
      END SUBROUTINE NEWDIA
