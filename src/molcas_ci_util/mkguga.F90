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
      SUBROUTINE MKGUGA(STSYM,Skip_MKSGNUM)
!
!     PURPOSE: MAKE THE GUGA TABLES
!     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!              THE START ADRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!              THREE USER DEFINED TYPES. Consult the struct.F90 and gugx.F90 files for the details.
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: IFCAS, CIS, SGS, EXS

      IMPLICIT None
!
      Integer STSYM
     Logical, Optional:: Skip_MKSGNUM

      Integer, Allocatable:: V11(:)
      Integer IAC, NVERT0

!Note that we do not associate the arrays here since the are not allocated yet.
      Associate (nVert => SGS%nVert, MidLev => SGS%MidLev, MVSta => SGS%MVSta,  &
     &           MVEnd => SGS%MVEnd,  nMidV => CIS%nMidV, nIpWlk => CIS%nIpWlk, &
     &           nWalk => CIS%nWalk,  MxUp => SGS%MxUp, MxDwn => SGS%MxDwn, &
     &           LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,               &
     &           LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,               &
     &           IA0 => SGS%IA0, IB0 => SGS%IB0, IC0 => SGS%IC0,       &
     &           nLev=>SGS%nLev, nSym=>SGS%nSym)
!
!     SET UP A FULL PALDUS DRT TABLE:
!     (INITIALLY NO RESTRICTIONS ARE PUT UP)
!
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6

      IF(IFCAS.NE.0) THEN
         CALL mma_allocate(SGS%DRT0,NVERT0,5,Label='DRT0')
         CALL mma_allocate(SGS%DOWN0,[1,NVERT0],[0,3],Label='DOWN0')
         SGS%DRTP => SGS%DRT0
         SGS%DOWNP=> SGS%DOWN0
      ELSE
         NVERT=NVERT0
         CALL mma_allocate(SGS%DRT,NVERT,5,Label='SGS%DRT')
         CALL mma_allocate(SGS%DOWN,[1,NVERT],[0,3],Label='SGS%DOWN')
         SGS%DRTP => SGS%DRT
         SGS%DOWNP=> SGS%DOWN
      ENDIF

      CALL mkDRT0(SGS,IA0,IB0,IC0,NVERT0,SGS%DRTP,SGS%DOWNP)
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' PALDUS DRT TABLE (UNRESTRICTED):'
      CALL PRDRT(SGS)
#endif
!
!     IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
!     VERTICES WHICH VIOLATE THE FORMER.
!
      IF(IFCAS.NE.0) THEN
        CALL mma_allocate(V11,NVERT0,Label='V11')
        CALL RESTR(NVERT0,SGS%DRT0,SGS%DOWN0,V11,LV1RAS, LV3RAS, LM1RAS, LM3RAS, NVERT)
!
!     REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
!
        CALL mma_allocate(SGS%DRT,nVert,5,Label='DRT')
        CALL mma_allocate(SGS%DOWN,[1,nVert],[0,3],Label='SGS%DOWN')
        CALL mkDRT(NVERT0,NVERT,SGS%DRT0,SGS%DOWN0,V11,SGS%DRT,SGS%DOWN)
        CALL mma_deallocate(V11)
        CALL mma_deallocate(SGS%DRT0)
        CALL mma_deallocate(SGS%DOWN0)
!
#ifdef _DEBUGPRINT_
        SGS%DOWNP=> SGS%DOWN
        SGS%DRTP => SGS%DRT
        Write(LF,*)
        Write(LF,*)' PALDUS DRT TABLE (RESTRICTED):'
        CALL PRDRT(SGS)
#endif

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
      CALL MKMID(SGS)
      NMIDV=MVEnd-MVSta+1
!
!     FORM VARIOUS OFFSET TABLES:
!     NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!           TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.
!
      CALL MKCOT(SGS,CIS)
!
!     CONSTRUCT THE CASE LIST
!
      Call MKCLIST(SGS,CIS)
!
!     SET UP ENUMERATION TABLES
!
      If (Present(Skip_MKSGNUM)) Then
         If (Skip_MKSGNUM) Return
      End If

      Call MKSGNUM(STSYM,SGS,CIS,EXS)

      End Associate
!
!     EXIT
!
      END SUBROUTINE MKGUGA

      SUBROUTINE MKGUGA_FREE()
!
!     PURPOSE: FREE THE GUGA TABLES
!
      use stdalloc, only: mma_deallocate
      use gugx, only:  SGS, CIS, EXS
      IMPLICIT None
!
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

      If (Allocated(CIS%NOW)) Call mma_deallocate(CIS%NOW)
      If (Allocated(CIS%IOW)) Call mma_deallocate(CIS%IOW)
      If (Allocated(CIS%NCSF)) Call mma_deallocate(CIS%NCSF)
      If (Allocated(CIS%NOCSF)) Call mma_deallocate(CIS%NOCSF)
      If (Allocated(CIS%IOCSF)) Call mma_deallocate(CIS%IOCSF)
      If (Allocated(CIS%ICASE)) Call mma_deallocate(CIS%ICASE)

      If (Allocated(EXS%NOCP)) Call mma_deallocate(EXS%NOCP)
      If (Allocated(EXS%IOCP)) Call mma_deallocate(EXS%IOCP)
      If (Allocated(EXS%ICOUP)) Call mma_deallocate(EXS%ICOUP)
      If (Allocated(EXS%VTAB)) Call mma_deallocate(EXS%VTAB)
      If (Allocated(EXS%MVL)) Call mma_deallocate(EXS%MVL)
      If (Allocated(EXS%MVR)) Call mma_deallocate(EXS%MVR)
      If (Allocated(EXS%SGTMP)) Call mma_deallocate(EXS%SGTMP)
      If (Allocated(EXS%USGN)) Call mma_deallocate(EXS%USGN)
      If (Allocated(EXS%LSGN)) Call mma_deallocate(EXS%LSGN)

      END SUBROUTINE MKGUGA_FREE
