************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,
     &                   nRs1,nRs2,nRs3,SGS,CIS,EXS,CIL,imode,ksym,
     &                   State_Sym)
*
      use stdalloc, only: mma_allocate, mma_deallocate
      use Str_Info, only: CFTP, CNSM
      use Struct, only: SGStruct, CIStruct, EXStruct
      Implicit None
      Integer nSym,iSpin,nActEl,nHole1,nElec3
      Integer nRs1(nSym), nRs2(nSym), nRs3(nSym)
      Type(SGStruct) SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS
      Real*8 CIL(*)
      Integer imode, ksym, State_Sym
*
      Real*8, Allocatable:: CINEW(:)
#ifdef _DEBUGPRINT_
      Real*8 :: PRWTHR=0.05d0
      Integer NICASE
#endif
      Integer nRas1T, nRas2T, nRas3T, iss, iS
      Integer NCONF
*
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

*
      Associate ( nLev=> SGS%nLev, nMidV =>CIS%nMidV,
     &            nVert =>SGS%nVert, MidLev=>SGS%MidLev,
     &            MVSta =>SGS%MVSta, MVEnd =>SGS%MVEnd,
     &            nIpWlk=>CIS%nIpWlk,
     &            MxUp => SGS%MxUp, MxDwn => SGS%MxDwn,
     &            LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,
     &            LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,
     &            IFRAS=>SGS%IFRAS)

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl
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

      Call mkGUGA(SGS,CIS)

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
      Call MKSGNUM(kSYM,SGS,CIS,EXS)

      nConf = CIS%nCSF(kSym)

      iss=1
      if (ksym.ne.state_sym) iss=2
*
#ifdef _DEBUGPRINT_
      NICASE = SIZE(CIS%ICASE)
      WRITE(6,101)
101   FORMAT(/,6X,100('-'),/,
     &      6X,29X,'Wave function printout: Split Graph format',/,
     &      6X, 8X,'in paranthesis: midvertex, upper-walk symmetry',
     &             ' upper- and lower-walk serial numbers',/,
     &         6X,100('-'),/)
      WRITE(6,102) PRWTHR
102   FORMAT(6X,'printout of CI-coefficients larger than',F6.2)
      Call SGPRWF_MCLR(ksym,PRWTHR,nSym,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,
     &                 NICASE,SGS%ISM,CIS%NOCSF,CIS%IOCSF,CIS%NOW,
     &                 CIS%IOW,CIS%ICASE,CIL)
      WRITE(6,103)
103   FORMAT(/,6X,100('-'),/)
#endif
*
      Call mma_allocate(CInew,NCONF,Label='CINew')
      Call REORD(NLEV,NVERT,MIDLEV,MVSta,NMIDV,MXUP,MXDWN,
     &           SGS%DRT,SGS%DOWN,SGS%DAW,SGS%UP,SGS%RAW,
     &           EXS%USGN,EXS%LSGN,nActEl,NLEV,NCONF,
     &           iMode,CNSM(iss)%ICONF,CFTP,kSym,
     &           CIL,CInew)
      Call DCopy_(nConf,CINew,1,CIL,1)
      Call mma_deallocate(CINew)

#ifdef _DEBUGPRINT_
      Call SGPRWF_MCLR(ksym,PRWTHR,nSym,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,
     &                 NICASE,SGS%ISM,CIS%NOCSF,CIS%IOCSF,
     &                 CIS%NOW,CIS%IOW,CIS%ICASE,CIL)
#endif
*
      End Associate

      End Subroutine GugaNew
