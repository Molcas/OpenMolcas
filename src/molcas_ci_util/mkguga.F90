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

contains

Subroutine SGClose(SGS)
use stdalloc, only: mma_deallocate
use Struct, only: SGStruct
Implicit None
Type (SGStruct) SGS

If (Allocated(SGS%ISM)) Call mma_deallocate(SGS%ISM)
If (Allocated(SGS%DRT0)) Call mma_deallocate(SGS%DRT0)
If (Allocated(SGS%DOWN0)) Call mma_deallocate(SGS%DOWN0)
If (Allocated(SGS%DRT)) Call mma_deallocate(SGS%DRT)
If (Allocated(SGS%DOWN)) Call mma_deallocate(SGS%DOWN)
If (Allocated(SGS%UP)) Call mma_deallocate(SGS%UP)
If (Allocated(SGS%MAW)) Call mma_deallocate(SGS%MAW)
If (Allocated(SGS%LTV)) Call mma_deallocate(SGS%LTV)
If (Allocated(SGS%DAW)) Call mma_deallocate(SGS%DAW)
If (Allocated(SGS%RAW)) Call mma_deallocate(SGS%RAW)
If (Allocated(SGS%SCR)) Call mma_deallocate(SGS%SCR)
If (Allocated(SGS%Ver)) Call mma_deallocate(SGS%Ver)
SGS%DRTP => Null()
SGS%DOWNP => Null()

end Subroutine SGClose

Subroutine CXClose(CIS,EXS)
use stdalloc, only: mma_deallocate
use Struct, only: CIStruct, EXStruct
Type (CIStruct) CIS
Type (ExStruct) ExS

IF (Allocated(CIS%NOW)) Call mma_deallocate(CIS%NOW)
IF (Allocated(CIS%IOW)) Call mma_deallocate(CIS%IOW)
IF (Allocated(CIS%NCSF)) Call mma_deallocate(CIS%NCSF)
IF (Allocated(CIS%NOCSF)) Call mma_deallocate(CIS%NOCSF)
IF (Allocated(CIS%IOCSF)) Call mma_deallocate(CIS%IOCSF)
IF (Allocated(CIS%ICase)) Call mma_deallocate(CIS%ICase)

IF (Allocated(EXS%NOCP)) Call mma_deallocate(EXS%NOCP)
IF (Allocated(EXS%IOCP)) Call mma_deallocate(EXS%IOCP)
IF (Allocated(EXS%ICoup)) Call mma_deallocate(EXS%ICoup)
IF (Allocated(EXS%VTab)) Call mma_deallocate(EXS%VTab)
IF (Allocated(EXS%SGTMP)) Call mma_deallocate(EXS%SGTMP)
IF (Allocated(EXS%MVL)) Call mma_deallocate(EXS%MVL)
IF (Allocated(EXS%MVR)) Call mma_deallocate(EXS%MVR)
If (Allocated(EXS%USGN)) Call mma_deallocate(EXS%USGN)
If (Allocated(EXS%LSGN)) Call mma_deallocate(EXS%LSGN)


end subroutine CXClose

      END SUBROUTINE MKGUGA_FREE
