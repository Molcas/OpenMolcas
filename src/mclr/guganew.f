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
      Subroutine GugaNew(ipCIL,imode,ksym)
*
      Implicit Real*8 (A-H,O-Z)
      Integer A0,B0,C0
*

#include "Input.fh"
#include "WrkSpc.fh"
#include "detdim.fh"
#include "csfbas_mclr.fh"
#include "spinfo_mclr.fh"
      Integer OrbSym(2*mxBas)
      Parameter (iPrint=0)
*
*     Call qEnter('GugaCtl')
*
      PRWTHR=0.05d0
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
      IAC=MIN(A0,C0)
      NVERT0=((A0+1)*(C0+1)*(2*B0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
      NDRT0=5*NVERT0
      NDOWN0=4*NVERT0
      NTMP=((NLEV+1)*(NLEV+2))/2
      Call GetMem('DRT0','ALLO','INTEGER',LDRT0,NDRT0)
      Call GetMem('DOWN','ALLO','INTEGER',LDOWN0,NDOWN0)
      Call GetMem('LTMP','ALLO','INTEGER',LTMP,NTMP)
      Call DRT0
     &     (A0,B0,C0,NVERT0,iWork(LDRT0),iWork(LDOWN0),
     &      NTMP,IWORK(LTMP))
      If ( iPrint.ge.5 ) Call PRDRT(NVERT0,iWork(LDRT0),iWork(LDOWN0))
      Call GetMem('LTMP','FREE','INTEGER',LTMP,NTMP)
*
      LV1RAS=ntRas1
      LV3RAS=LV1RAS+ntRas2
      LM1RAS=2*LV1RAS-nHole1
      LM3RAS=nActEl-nElec3
      Call GetMem('LV11','ALLO','INTEGER',LV,NVERT0)
      Call RESTR
     &     (NVERT0,iWork(LDRT0),iWork(LDOWN0),iWork(LV),
     &      LV1RAS,LV3RAS,LM1RAS,LM3RAS,NVERT)
*
      NDRT=5*NVERT
      NDOWN=4*NVERT
      Call GetMem('DRT1','ALLO','INTEGER',LDRT,NDRT)
      Call GetMem('DWN1','ALLO','INTEGER',LDOWN,NDOWN)
      Call DRT
     &     (NVERT0,NVERT,iWork(LDRT0),iWork(LDOWN0),iWork(LV),
     &      iWork(LDRT),iWork(LDOWN))
      If ( iPrint.ge.5 ) Call PRDRT(NVERT,iWork(LDRT),iWork(LDOWN))
      Call GetMem('LV11','FREE','INTEGER',LV,NVERT0)
      Call GetMem('DRT0','FREE','INTEGER',LDRT0,NDRT0)
      Call GetMem('DOWN','FREE','INTEGER',LDOWN0,NDOWN0)
*
      NDAW=5*NVERT
      Call GetMem('DAW1','ALLO','INTEGER',LDAW,NDAW)
      Call MKDAW(NVERT,iWork(LDOWN),iWork(LDAW),iPrint)
*
      NUP=4*NVERT
      NRAW=5*NVERT
      Call GetMem('LUP1','ALLO','INTEGER',LUP,NUP)
      Call GetMem('RAW1','ALLO','INTEGER',LRAW,NRAW)
      Call MKRAW
     &     (NVERT,iWork(LDOWN),iWork(LDAW),iWork(LUP),iWork(LRAW),
     &           iPrint)
*
      NLTV=NLEV+2
      Call GetMem('LTV1','ALLO','INTEGER',LLTV,NLTV)
      Call MKMID
     &     (NVERT,NLEV,iWork(LDRT),
     &      iWork(LDOWN),iWork(LDAW),iWork(LUP),iWork(LRAW),
     &      iWork(LLTV),
     &      MIDLEV,NMIDV,MIDV1,MIDV2,MXUP,MXDWN,iPrint)
      Call GetMem('LTV1','FREE','INTEGER',LLTV,NLTV)
*
      NIPWLK=1+(MIDLEV-1)/15
      NIPWLK=MAX(NIPWLK,1+(NLEV-MIDLEV-1)/15)
      NNOW=2*NMIDV*nSym
      NIOW=NNOW
      NNOCSF=NMIDV*(nSym**2)
      NIOCSF=NNOCSF
      NSCR=3*(NLEV+1)
      Call GetMem('NOW1','ALLO','INTEGER',LNOW,NNOW)
      Call GetMem('IOW1','ALLO','INTEGER',LIOW,NIOW)
      Call GetMem('NCSF','ALLO','INTEGER',LNOCSF,NNOCSF)
      Call GetMem('ICSF','ALLO','INTEGER',LIOCSF,NIOCSF)
      Call GetMem('SCR1','ALLO','INTEGER',LSCR,NSCR)
*     Call GetMem('NCSF','ALLO','INTEGER',LNCSF,nSym)
      Call MKCOT
     &     (nSym,NLEV,NVERT,MIDLEV,NMIDV,MIDV1,MIDV2,NWALK,NIPWLK,
     &      OrbSym,iWork(LDOWN),iWork(LNOW),iWork(LIOW),
     &      NCSF,iWork(LIOCSF),iWork(LNOCSF),
     &      iWork(LSCR),iPrint)
*
*     If ( nConf.ne.NCSF(state_sym) ) then
*        Write (*,*)
*        Write (*,*) ' *** Error in subroutine GUGACTL ***'
*        Write (*,*) ' Inconsistent, number of configurations'
*        Write (*,*)
*     End If
      NCONF=NCSF(kSym)
      iss=1
      if (ksym.ne.state_sym) iss=2
*
      NILNDW=NWALK
      NICASE=NWALK*NIPWLK
      Call GetMem('CASE','ALLO','INTEG',LICASE,NICASE)
      Call MKCLIST_MCLR
     &     (nSym,NLEV,NVERT,MIDLEV,MIDV1,MIDV2,NMIDV,NICASE,NIPWLK,
     &      OrbSym,iWork(LDOWN),iWork(LNOW),iWork(LIOW),
     &      iWork(LICASE),iWork(LSCR),iPrint)
      Call GetMem('SCR1','FREE','INTEG',LSCR,NSCR)
*
      NUSGN=MXUP*NMIDV
      NLSGN=MXDWN*NMIDV
      Call GetMem('IUSG','ALLO','INTEG',LUSGN,NUSGN)
      Call GetMem('ILSG','ALLO','INTEG',LLSGN,NLSGN)
      Call MKSGNUM
     &     (ksym,nSym,NLEV,
     &       NVERT,MIDLEV,NMIDV,MXUP,MXDWN,NICASE,NIPWLK,
     &      iWork(LDOWN),iWork(LUP),iWork(LDAW),iWork(LRAW),
     &      iWork(LNOW),iWork(LIOW),iWork(LUSGN),iWork(LLSGN),
     &      iWork(LICASE),iPrint)
*
      If (iPrint.ge.5) Then
      If (imode.eq.0.and.iAnd(kprint,8).eq.8) Then
      WRITE(6,101)
101   FORMAT(/,6X,100(1H-),/,
     &      6X,29X,'Wave function printout: Split Graph format',/,
     &      6X, 8X,'in paranthesis: midvertex, upper-walk symmetry',
     &             ' upper- and lower-walk serial numbers',/,
     &         6X,100(1H-),/)
      WRITE(6,102) PRWTHR
102   FORMAT(6X,'printout of CI-coefficients larger than',F6.2)
      Call SGPRWF
     &     (ksym,PRWTHR,
     &      nSym,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,NICASE,
     &      OrbSym,iWork(LNOCSF),iWork(LIOCSF),
     &      iWork(LNOW),iWork(LIOW),
     &      iWork(LICASE),Work(ipCIL))
      WRITE(6,103)
103   FORMAT(/,6X,100(1H-),/)
      End If
      End If
*
      jPrint=iPrint
      If (TimeDep) Then
         Call GetMem('CIvec','Allo','Real',ipCInew,2*NCONF)
      Else
         Call GetMem('CIvec','Allo','Real',ipCInew,NCONF)
      End If
      Call REORD
     &     (NLEV,NVERT,MIDLEV,MIDV1,MIDV2,NMIDV,MXUP,MXDWN,
     &      iWork(LDRT),iWork(LDOWN),iWork(LDAW),
     &      iWork(LUP),iWork(LRAW),iWork(LUSGN),iWork(LLSGN),
     &      nActEl,NLEV,NCONF,NTYP,
     &      iMode,jPrint,
     &      iWork(KICONF(iss)),
     &      iWork(KCFTP),NCNATS(1,kSym),NCPCNT,
     &      Work(ipCIL),Work(ipCInew),minop)
      If (imode.eq.0.and.iAnd(kprint,8).eq.8)
     &Call SGPRWF
     &     (ksym,PRWTHR,
     &      nSym,NLEV,NCONF,MIDLEV,NMIDV,NIPWLK,NICASE,
     &      OrbSym,iWork(LNOCSF),iWork(LIOCSF),
     &      iWork(LNOW),iWork(LIOW),
     &      iWork(LICASE),Work(ipCInew))
CEAW970812     &      Work(ipCIL),Work(ipCInew))
      If (TimeDep) Then
         Call GetMem('OCIvec','Free','Real',ipCIL,2*nConf)
      Else
         Call GetMem('OCIvec','Free','Real',ipCIL,nConf)
      End If
      ipCIL = ipCInew
*
      Call GetMem('ILSG','FREE','INTEGER',LLSGN,NLSGN)
      Call GetMem('IUSG','FREE','INTEGER',LUSGN,NUSGN)
      Call GetMem('CASE','FREE','INTEGER',LICASE,NICASE)
*     Call GetMem('NCSF','FREE','INTEGER',LNCSF,nSym)
      Call GetMem('ICSF','FREE','INTEGER',LIOCSF,NIOCSF)
      Call GetMem('NCSF','FREE','INTEGER',LNOCSF,NNOCSF)
      Call GetMem('IOW1','FREE','INTEGER',LIOW,NIOW)
      Call GetMem('NOW1','FREE','INTEGER',LNOW,NNOW)
      Call GetMem('RAW1','FREE','INTEGER',LRAW,NRAW)
      Call GetMem('LUP1','FREE','INTEGER',LUP,NUP)
      Call GetMem('DAW1','FREE','INTEGER',LDAW,NDAW)
      Call GetMem('DWN1','FREE','INTEGER',LDOWN,NDOWN)
      Call GetMem('DRT1','FREE','INTEGER',LDRT,NDRT)
*
*     Call qExit('GugaCtl')
*
      Return
      End
