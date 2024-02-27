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
      use gugx, only: IFCAS, LV1RAS, LM1RAS, LV3RAS, LM3RAS, CIS, SGS

      IMPLICIT REAL*8 (A-H,O-Z)
!
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
      Character(LEN=16), Parameter :: ROUTINE='GUGACTL '
      Integer nVert0

      Interface
      SUBROUTINE MKGUGA(NLEV,NSYM,STSYM,Skip_MKSGNUM)
      IMPLICIT None

      Integer NLEV, NSYM, STSYM
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

      Associate( nLev =>SGS%nLev,   &
     &            IA0 => SGS%IA0, IB0 => SGS%IB0, IC0 => SGS%IC0)

! Local print level (if any)
#ifdef _DEBUGPRINT_
      WRITE(LF,*)' Entering ',ROUTINE
#endif
!
!     SET IFCAS FLAG
!     IFCAS = 0 : THIS IS A CAS CALCULATION
!     IFCAS = 1 : THIS IS A RAS CALCULATION
!
      IFCAS=0
      IF (NHOLE1.NE.0.OR.NELEC3.NE.0) IFCAS=1
      DO IS=1,NSYM
        IF (IFCAS.NE.0.AND.NASH(IS).NE.0)IFCAS=IFCAS+1
      END DO
!
!     CREATE THE SYMMETRY INDEX VECTOR
!
      CALL MKNSM()
!
!     (IFCAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
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

      Call mma_allocate(SGS%ISM,nLev,Label='SGS%ISM')

      SGS%ISM(1:nLev)=NSM(1:nLev)
!
!     COMPUTE RAS RESTRICTIONS ON VERTICES:
!
      LV3RAS=LV1RAS+NRAS2
      LM1RAS=2*LV1RAS-NHOLE1
      LM3RAS=NACTEL-NELEC3
!
!     COMPUTE TOP ROW OF THE GUGA TABLE
!
      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=NLEV-IA0-IB0
      IF ( ((2*IA0+IB0).NE.NACTEL) .OR.                                 &
     &     (IA0.LT.0) .OR.                                              &
     &     (IB0.LT.0) .OR.                                              &
     &     (IC0.LT.0) ) then
        Write(LF,*)'GUGACTL Error: Impossible specifications.'
        Write(LF,'(1x,a,3I8)')'NACTEL,NLEV,ISPIN:',NACTEL,NLEV,ISPIN
        Write(LF,'(1x,a,3I8)')'IA0,IB0,IC0:      ',IA0,IB0,IC0
        Write(LF,*)' This is a severe internal error, or possibly'
        Write(LF,*)' indicates a strange input which should have been'
        Write(LF,*)' diagnosed earlier. Please submit a bug report.'
        Call Quit(_RC_GENERAL_ERROR_)
      End If
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
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
      CALL MKGUGA(NLEV,NSYM,STSYM)
      NCONF=CIS%NCSF(STSYM)
      If ( NAC.eq.0 ) NCONF=1

      End Associate

      END SUBROUTINE GUGACTL
