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
! Copyright (C) Per Ake Malmqvist                                      *
!               Markus P. Fuelscher                                    *
!***********************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE GUGACTL()
!
!     PURPOSE: CONTROL ROUTINE TO SET UP GUGA TABLES
!     AUTHOR:  P.-AA. MALMQVIST
!
!     MODIFIED TO FIT THE DETRAS PROGRAM BY M.P. FUELSCHER
!
      use Definitions, only: LF => u6
      use stdalloc, only: mma_allocate
      use gugx, only: IFRAS, CIS, SGS

      IMPLICIT REAL*8 (A-H,O-Z)
!
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
      Character(LEN=16), Parameter :: ROUTINE='GUGACTL '

      Interface
      SUBROUTINE MKGUGA(STSYM,Skip_MKSGNUM)
      IMPLICIT None

      Integer STSYM
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl

      Associate( nLev =>SGS%nLev,   &
     &           LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,   &
     &           LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,   &
     &           IA0 => SGS%IA0, IB0 => SGS%IB0, IC0 => SGS%IC0, &
     &           nVert0 => SGS%nVert0)

! Local print level (if any)
#ifdef _DEBUGPRINT_
      WRITE(LF,*)' Entering ',ROUTINE
#endif
!
!     SET IFRAS FLAG
!     IFRAS = 0 : THIS IS A CAS CALCULATION
!     IFRAS = 1 : THIS IS A RAS CALCULATION
!
      IFRAS=0
      IF (NHOLE1.NE.0.OR.NELEC3.NE.0) IFRAS=1
      DO IS=1,NSYM
        IF (IFRAS.NE.0.AND.NASH(IS).NE.0)IFRAS=IFRAS+1
      END DO
!
!     CREATE THE SYMMETRY INDEX VECTOR
!
      CALL MKISM(SGS)
!
!     (IFRAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
!     IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
!     SET UP THE REINDEXING TABLE
!
      CALL SETSXCI()
!
!     FIND TOTAL NUMBER OF VERTICES IN THE SUBSPACES
!
!... for RAS
      NLEV=0
      DO IS=1,NSYM
        NLEV=NLEV+NRS1(IS)
      END DO
      LV1RAS=NLEV
      DO IS=1,NSYM
        NLEV=NLEV+NRS2(IS)
      END DO
      NRAS2=NLEV-LV1RAS
      DO IS=1,NSYM
        NLEV=NLEV+NRS3(IS)
      END DO

!
!     COMPUTE RAS RESTRICTIONS ON VERTICES:
!
      LV3RAS=LV1RAS+NRAS2
      LM1RAS=2*LV1RAS-NHOLE1
      LM3RAS=NACTEL-NELEC3
!
!     COMPUTE TOP ROW OF THE GUGA TABLE
!
      Call mknVert0(SGS)

      If ( NVERT0.eq.0 ) then
        NCONF=0
        Return
      End If
      If ( doBlockDMRG ) then
        NCONF=1
        Return
      End If
!
!     INITIALIZE GUGA TABLES:
!
      CALL MKGUGA(STSYM)

      NCONF=CIS%NCSF(STSYM)
      If ( NAC.eq.0 ) NCONF=1

      End Associate

      END SUBROUTINE GUGACTL
