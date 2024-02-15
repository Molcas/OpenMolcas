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
      Subroutine GugaNew(CIL,imode,ksym)
*
      use gugx, only: NLEV, A0 => IA0, B0 => IB0, C0 => IC0,
     &                NVERT,MIDLEV,MIDV1,MIDV2,NMIDV,MXUP,MXDWN,DRT,
     &                DOWN,DAW,UP,RAW,USGN,LSGN,ICASE,IFCAS,
     &                LV1RAS, LV3RAS, LM1RAS, LM3RAS, NOCSF, IOCSF,
     &                NOW => NOW1, IOW => IOW1, NICASE, NIPWLK
      use Str_Info, only: CFTP, CNSM
      Implicit None
      Integer imode, ksym
      Real*8 CIL(*)
*
#include "Input.fh"
#include "stdalloc.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
      Integer OrbSym(2*mxBas)
      Integer, Parameter:: iPrint=0
      Real*8, Allocatable:: CINEW(:)
      Real*8 :: PRWTHR=0.05d0
      Integer ntRas1, ntRas2, ntRas3, iSym, jPrint, iBas, iOrb, iss
*
*
      ntRas1=0
      ntRas2=0
      ntRas3=0
      Do iSym=1,nSym
         ntRas1=ntRas1+nRs1(iSym)
         ntRas2=ntRas2+nRs2(iSym)
         ntRas3=ntRas3+nRs3(iSym)
      End Do
*
      B0=iSpin-1
      A0=(nActEl-B0)/2
      C0=ntASh-A0-B0
      If ( (2*A0+B0).ne.nActEl ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' 2*A0+B0.ne.nActEl '
         Write (6,*)
      End If
      If ( A0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' A0.lt.0'
         Write (6,*)
      End If
      If ( B0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' B0.lt.0'
         Write (6,*)
      End If
      If ( C0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' C0.lt.0'
         Write (6,*)
      End If
*
      iOrb=0
      Do iSym=1,nSym
         Do iBas=1,nRs1(iSym)
            iOrb=iOrb+1
            OrbSym(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs2(iSym)
            iOrb=iOrb+1
            OrbSym(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs3(iSym)
            iOrb=iOrb+1
            OrbSym(iOrb)=iSym
         End Do
      End Do
*
      NLEV=ntASh
      LV1RAS=ntRas1
      LV3RAS=LV1RAS+ntRas2
      LM1RAS=2*LV1RAS-nHole1
      LM3RAS=nActEl-nElec3

      IFCAS=1
      Call mkGUGA(OrbSym,NSYM,kSym,NCSF)
      NCONF=NCSF(kSym)
      iss=1
      if (ksym.ne.state_sym) iss=2
*
      If (iPrint.ge.5) Then
      If (imode.eq.0.and.iAnd(kprint,8).eq.8) Then
      WRITE(6,101)
101   FORMAT(/,6X,100('-'),/,
     &      6X,29X,'Wave function printout: Split Graph format',/,
     &      6X, 8X,'in paranthesis: midvertex, upper-walk symmetry',
     &             ' upper- and lower-walk serial numbers',/,
     &         6X,100('-'),/)
      WRITE(6,102) PRWTHR
102   FORMAT(6X,'printout of CI-coefficients larger than',F6.2)
      Call SGPRWF_MCLR(ksym,PRWTHR,nSym,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,
     &                 NICASE,OrbSym,NOCSF,IOCSF,NOW,IOW,ICASE,CIL)
      WRITE(6,103)
103   FORMAT(/,6X,100('-'),/)
      End If
      End If
*
      jPrint=iPrint
      Call mma_allocate(CInew,NCONF,Label='CINew')
      Call REORD(NLEV,NVERT,MIDLEV,MIDV1,MIDV2,NMIDV,MXUP,MXDWN,DRT,
     &           DOWN,DAW,UP,RAW,USGN,LSGN,nActEl,NLEV,NCONF,NTYP,iMode,
     &           jPrint,CNSM(iss)%ICONF,CFTP,NCNATS(1,kSym),NCPCNT,
     &           CIL,CInew,minop)
      If (imode.eq.0.and.iAnd(kprint,8).eq.8)
     &Call SGPRWF_MCLR(ksym,PRWTHR,nSym,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,
     &                 NICASE,OrbSym,NOCSF,IOCSF,NOW,IOW,ICASE,CInew)
      Call DCopy_(nConf,CINew,1,CIL,1)
      Call mma_deallocate(CINew)
*
      Call mkGUGA_Free()

      End Subroutine GugaNew
