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
      Subroutine GugaNew(SGS,CIS,EXS,CIL,imode,ksym)
*
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: IFRAS
      use Str_Info, only: CFTP, CNSM
      use Struct, only: SGStruct, CIStruct, EXStruct
      Implicit None
      Type(SGStruct) SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS
      Integer imode, ksym
      Real*8 CIL(*)
*
#include "Input.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
      Real*8, Allocatable:: CINEW(:)
      Real*8 :: PRWTHR=0.05d0
      Integer ntRas1, ntRas2, ntRas3, iSym, iss
      Integer NICASE
*
      Interface
      SUBROUTINE MKGUGA(STSYM,Skip_MKSGNUM)
      IMPLICIT None

      Integer STSYM
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

*
      NICASE = SIZE(CIS%ICASE)

      Associate ( nLev=> SGS%nLev, nMidV =>CIS%nMidV,
     &            nVert =>SGS%nVert, MidLev=>SGS%MidLev,
     &            MVSta =>SGS%MVSta, MVEnd =>SGS%MVEnd,
     &            nIpWlk=>CIS%nIpWlk,
     &            MxUp => SGS%MxUp, MxDwn => SGS%MxDwn,
     &            LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,
     &            LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS)

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl
*
      Call mkISm(SGS)

      Call mknVert0(SGS)
*
      ntRas1=0
      ntRas2=0
      ntRas3=0
      Do iSym=1,nSym
         ntRas1=ntRas1+nRs1(iSym)
         ntRas2=ntRas2+nRs2(iSym)
         ntRas3=ntRas3+nRs3(iSym)
      End Do

      LV1RAS=ntRas1
      LV3RAS=LV1RAS+ntRas2
      LM1RAS=2*LV1RAS-nHole1
      LM3RAS=nActEl-nElec3

      IFRAS=1

      Call mkGUGA(kSym)

      NCSF(1:nSym)=CIS%NCSF(1:nSym)
      NCONF=CIS%NCSF(kSym)

      iss=1
      if (ksym.ne.state_sym) iss=2
*
#ifdef _DEBUGPRINT_
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
     &           EXS%USGN,EXS%LSGN,
     &           nActEl,NLEV,NCONF,NTYP,
     &           iMode,CNSM(iss)%ICONF,CFTP,NCNATS(1,kSym),
     &           NCPCNT,CIL,CInew,minop)
      If (imode.eq.0.and.iAnd(kprint,8).eq.8)
     &Call SGPRWF_MCLR(ksym,PRWTHR,nSym,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,
     &                 NICASE,SGS%ISM,CIS%NOCSF,CIS%IOCSF,
     &                 CIS%NOW,CIS%IOW,CIS%ICASE,CInew)
      Call DCopy_(nConf,CINew,1,CIL,1)
      Call mma_deallocate(CINew)
*
      Call mkGUGA_Free()

      End Associate

      End Subroutine GugaNew
