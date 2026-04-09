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
* Copyright (C) 1989,1998, Per Ake Malmqvist                           *
************************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE CITRA     IBM-3090 RELEASE 89 01 31
*  USE THE COEFFICIENTS FOR A SEQUENCE OF SINGLE-ORBITAL TRANSFOR-
*  MATION, TRA, TO TRANSFORM THE CI EXPANSION COEFFICIENTS
*  IN-PLACE TO A NEW NON-ON ORBITAL BASIS.
*  NEW VERSION 981122, using user define types SGS,CIS,XS.
************************************************************************
*  CITRA
*
*> @brief
*>   Recompute a CI coefficient array to use another orbital basis
*> @author P. &Aring;. Malmqvist
*>
*> @details
*> For a given linear transformation of the orbitals, and a
*> CI array where the CSF basis is built from the old
*> orbitals, compute the CI array using the new orbitals
*> instead. The orbitals are transformed sequentially, and
*> for each active orbital, a call to ::SSOTRA performs the
*> single-orbital transformation.
*>
*> @param[in]     WFTP Wave function Type Name
*> @param[in]     SGS Split Graph Structure user defined type
*> @param[in]     CIS CI Structure user define type
*> @param[in]     EXS  Excitation operator Structure user defined type
*> @param[in]     LSM  Wave function Symmetry Label
*> @param[in]     TRA  Transformation Matrix
*> @param[in]     NCO  Number of Configuration Functions
*> @param[in,out] CI   CI Array
************************************************************************
      SUBROUTINE CITRA(WFTP,SGS,CIS,EXS,LSM,TRA,NCO,CI)
      use definitions, only: iwp, wp
#ifdef DEBUG_MPSSI
      use definitions, only: u6
#endif
      use constants, only: One
      use gugx, only: SGStruct, CIStruct, EXStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NTRA,NOSH,NISH,NASH
      IMPLICIT NONE
      CHARACTER(LEN=8), Intent(in):: WFTP
      Type (SGStruct), intent(in):: SGS
      Type (CIStruct), intent(in):: CIS
      Type (EXStruct), intent(in):: EXS
      integer(kind=iwp), intent(in):: LSM, NCO
      Real(kind=wp), intent(in)::TRA(NTRA)
      Real(kind=wp), intent(inout)::CI(NCO)

      real(kind=wp), Allocatable:: TMP(:)
      real(kind=wp) FAC,CKK
      integer(kind=iwp) ISTA,ISYM,NO,I,II,NA,NI


#ifdef DEBUG_MPSSI
      write(6,*)' Entering CITRA. norm=',ddot_(NCO,CI,1,CI,1)
#endif
!     write(6,*)' Entering CITRA. TRA='
!     write(6,'(1x,5f16.8)')(TRA(I),I=1,NTRA)
!     write(6,*)' Entering CITRA. CI='
!     write(6,'(1x,5f16.8)')(CI(I),I=1,NCO)
C TRA contains square matrices, one per symmetry
C  FIRST TRANSFORM THE INACTIVE ORBITALS:
      FAC=One
      ISTA=1
      DO ISYM=1,NSYM
        NO=NOSH(ISYM)
        DO I=1,NISH(ISYM)
          II=ISTA+(NO+1)*(I-1)
          CKK=TRA(II)
          FAC=FAC*CKK
        END DO
        ISTA=ISTA+NO**2
      END DO
!     write(6,*) 'FAC, FAC**2 ... ',FAC,FAC**2
      FAC=FAC**2
      CALL DSCAL_(NCO,FAC,CI,1)
!     write(6,*)' CITRA. inactive done CI='
!     write(6,'(1x,5f16.8)')(CI(I),I=1,NCO)
C  THEN THE ACTIVE ONES:
      if (WFTP /= 'EMPTY   ') then
* The HISPIN case may be buggy and is not presently used.
      IF(WFTP.EQ.'HISPIN  '.or.WFTP.EQ.'CLOSED  ') THEN
        ISTA=1
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NO=NOSH(ISYM)
          DO I=NI+1,NI+NA
            II=ISTA+(NO+1)*(I-1)
            CKK=TRA(II)
            FAC=FAC*CKK
          END DO
          ISTA=ISTA+NO**2
        END DO
        IF(WFTP.EQ.'CLOSED  ') FAC=FAC**2
        CALL DSCAL_(NCO,FAC,CI,1)
      ELSE
C The general case:
        CALL mma_allocate(TMP,NCO,Label='TMP')
        ISTA=1
        DO ISYM=1,NSYM
          NA=NASH(ISYM)
          NO=NOSH(ISYM)
          IF(NA.NE.0) THEN
            CALL SSOTRA(SGS,CIS,EXS,ISYM,LSM,NA,NO,
     &                TRA(ISTA),NCO,CI,TMP)
          END IF
          ISTA=ISTA+NO**2
        END DO
        CALL mma_deallocate(TMP)
      END IF
#ifdef DEBUG_MPSSI
      write(6,*)' DONE in  CITRA. norm=',ddot_(NCO,CI,1,CI,1)
#endif
!     write(6,*)' CITRA completely done. CI='
!     write(6,'(1x,5f16.8)')(CI(I),I=1,NCO)

      end if

      END SUBROUTINE CITRA
