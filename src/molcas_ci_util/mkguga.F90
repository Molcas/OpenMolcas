!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE MKGUGA(SGS,CIS)
!
!     PURPOSE: MAKE THE GUGA TABLES
!     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!              THE START ADRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!              THREE USER DEFINED TYPES. Consult the struct.F90 and gugx.F90 files for the details.
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate
      use struct, only: SGStruct, CIStruct
      IMPLICIT None
      Type(SGStruct),Target:: SGS
      Type(CIStruct) CIS

      Interface
      Subroutine MKDRT0(SGS)
      use struct, only: SGStruct
      IMPLICIT None
      Type(SGStruct), Target:: SGS
      End Subroutine MKDRT0
      End Interface
!
!     CREATE THE SYMMETRY INDEX VECTOR
!
      CALL MKISM(SGS)
!
!     COMPUTE TOP ROW OF THE GUGA TABLE
!
      Call mknVert0(SGS)

!Note that we do not associate the arrays here since the are not allocated yet.
      Associate (nVert => SGS%nVert, nVert0=>SGS%nVert0, IFRAS=>SGS%IFRAS)
!
!     SET UP A FULL PALDUS DRT TABLE:
!     (INITIALLY NO RESTRICTIONS ARE PUT UP)
!
      NVERT=NVERT0
      IF(IFRAS.NE.0) THEN
         CALL mma_allocate(SGS%DRT0,NVERT0,5,Label='DRT0')
         CALL mma_allocate(SGS%DOWN0,[1,NVERT0],[0,3],Label='DOWN0')
         SGS%DRTP => SGS%DRT0
         SGS%DOWNP=> SGS%DOWN0
      ELSE
         CALL mma_allocate(SGS%DRT,NVERT,5,Label='SGS%DRT')
         CALL mma_allocate(SGS%DOWN,[1,NVERT],[0,3],Label='SGS%DOWN')
         SGS%DRTP => SGS%DRT
         SGS%DOWNP=> SGS%DOWN
      ENDIF

      CALL mkDRT0(SGS)
!
!     IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
!     VERTICES WHICH VIOLATE THE FORMER.
!
      IF(IFRAS.NE.0) THEN
        CALL mkRAS(SGS)
!
!     REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
!
        CALL mkDRT(SGS)

!     IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED
!     DRT TABLE
!
      ENDIF

      SGS%DOWNP=> Null()
      SGS%DRTP => Null()
!
!     CALCULATE ARC WEIGHT.
!
      CALL MKDAW(SGS)
!
!     COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
!
      CALL MKRAW(SGS)
!
!     COMPUTE LTV TABLES.
!
      CALL MKLTV(SGS)
!
!     COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.
!
      CALL MKMID(SGS,CIS)
!
      End Associate
!
!     EXIT
!
      END SUBROUTINE MKGUGA

      SUBROUTINE MKGUGA_FREE(SGS,CIS,EXS)
!
!     PURPOSE: FREE THE GUGA TABLES
!
      use struct, only: SGStruct, CIStruct, EXStruct
      IMPLICIT None
      Type(SGStruct),Target:: SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS
!
      Call sgclose(SGS)

      Call cxclose(CIS,EXS)

      END SUBROUTINE MKGUGA_FREE
