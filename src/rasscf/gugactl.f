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
* Copyright (C) Per Ake Malmqvist                                      *
*               Markus P. Fuelscher                                    *
************************************************************************
      SUBROUTINE GUGACTL
C
C     PURPOSE: CONTROL ROUTINE TO SET UP GUGA TABLES
C     AUTHOR:  P.-AA. MALMQVIST
C
C     MODIFIED TO FIT THE DETRAS PROGRAM BY M.P. FUELSCHER
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "warnings.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='GUGACTL ')
#include "gugx.fh"
#include "WrkSpc.fh"

C Local print level (if any)
      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
C
C     SET IFCAS FLAG
C     IFCAS = 0 : THIS IS A CAS CALCULATION
C     IFCAS = 1 : THIS IS A RAS CALCULATION
C
      IFCAS=0
      IF (NHOLE1.NE.0.OR.NELEC3.NE.0) IFCAS=1
      DO IS=1,NSYM
        IF (IFCAS.NE.0.AND.NASH(IS).NE.0)IFCAS=IFCAS+1
      END DO
C
C     CREATE THE SYMMETRY INDEX VECTOR
C
      CALL MKNSM
C
C     (IFCAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
C     IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
C     SET UP THE REINDEXING TABLE
C
      CALL SETSXCI
C
C     FIND TOTAL NUMBER OF VERTICES IN THE SUBSPACES
C
C... for RAS
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
C
C     COMPUTE RAS RESTRICTIONS ON VERTICES:
C
      LV3RAS=LV1RAS+NRAS2
      LM1RAS=2*LV1RAS-NHOLE1
      LM3RAS=NACTEL-NELEC3
C
C     COMPUTE TOP ROW OF THE GUGA TABLE
C
      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=NLEV-IA0-IB0
      IF ( ((2*IA0+IB0).NE.NACTEL) .OR.
     &     (IA0.LT.0) .OR.
     &     (IB0.LT.0) .OR.
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
        Goto 100
      End If
      If ( doBlockDMRG ) then
        NCONF=1
        Goto 100
      End If
C
C     INITIALIZE GUGA TABLES:
C
      CALL MKGUGA(NSM,IPRLEV)
      NCONF=NCSF(STSYM)
      If ( NAC.eq.0 ) NCONF=1

100   Continue

      RETURN
      END
