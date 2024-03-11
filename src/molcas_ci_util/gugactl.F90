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
      SUBROUTINE GUGACTL(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,STSYM,DoblockDMRG)
!
!     PURPOSE: CONTROL ROUTINE TO SET UP GUGA TABLES
!     AUTHOR:  P.-AA. MALMQVIST
!
!     MODIFIED TO FIT THE DETRAS PROGRAM BY M.P. FUELSCHER
!
#ifdef _DEBUGPRINT_
      use Definitions, only: u6
#endif
      use Struct, only: SGStruct, CIStruct, EXStruct
      IMPLICIT None
      Integer nSym,iSpin,nActEl,nHole1,nElec3,STSYM
      Integer nRs1(nSym),nRs2(nSym),nRs3(nSym)
      Type(SGStruct) SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS
      Logical DoBlockDMRG
!
      Integer IS, nRas1T,nRas2T,nRas3T

      Interface
      SUBROUTINE MKGUGA(SGS,CIS)
      use Struct, only: SGStruct, CIStruct
      IMPLICIT None

      Type(SGStruct), Target:: SGS
      Type(CIStruct) CIS
      End SUBROUTINE MKGUGA
      End Interface

      nRas1T=Sum(nRs1(1:nSym))
      nRas2T=Sum(nRs2(1:nSym))
      nRas3T=Sum(nRs3(1:nSym))

#ifdef _DEBUGPRINT_
      Write (u6,*) 'nSym,iSpin,nActEl,nHole1,nElec3,nRas1T,nRas2T,nRas3T,STSYM=', &
                   nSym,iSpin,nActEl,nHole1,nElec3,nRas1T,nRas2T,nRas3T,STSYM
#endif

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl

      Associate( nLev =>SGS%nLev,   &
     &           LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,   &
     &           LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,   &
     &           nVert0 => SGS%nVert0, IFRAS=>SGS%IFRAS)

!
!     COMPUTE RAS RESTRICTIONS ON VERTICES:
!
      LV1RAS=NRAS1T
      LV3RAS=nRas1T+NRAS2T
      LM1RAS=2*nRas1T-NHOLE1
      LM3RAS=NACTEL-nElec3

!     SET IFRAS FLAG
!     IFRAS = 0 : THIS IS A CAS CALCULATION
!     IFRAS = 1 : THIS IS A RAS CALCULATION
!
      IF ((NRAS1T+NRAS3T)/=0) Then
         IFRAS=1
      Else
         IFRAS=0
      End If
      DO IS=1,NSYM
        IF (IFRAS.NE.0.AND.nRs2(IS).NE.0)IFRAS=IFRAS+1
      END DO
!
!     INITIALIZE GUGA TABLES:
!
      CALL MKGUGA(SGS,CIS)

      If ( NVERT0.eq.0 ) then
        CIS%NCSF(STSYM)=0
        Return
      End If
      If ( doBlockDMRG ) then
        CIS%NCSF(STSYM)=1
        Return
      End If

!     PURPOSE: FREE THE GUGA TABLES
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
      Call MKSGNUM(STSYM,SGS,CIS,EXS)

      If ( NActEl.eq.0 ) CIS%NCSF(STSYM)=1

      End Associate
!
!     (IFRAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
!     IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
!     SET UP THE REINDEXING TABLE
!
      CALL SETSXCI()

      END SUBROUTINE GUGACTL
