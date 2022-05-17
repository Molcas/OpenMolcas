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
      Subroutine GugaCtl_dmrg()
*

      Implicit Real*8 (A-H,O-Z)
      Integer A0,B0,C0
*
#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
      Integer OrbSym(2*mxBas)
      Integer, Parameter:: iPrint=0
      Integer, Allocatable:: DRT0(:), DOWN0(:), TMP(:), V11(:), DRT(:),
     &                       DOWN(:), DAW(:), UP(:), RAW(:), LTV(:),
     &                       NOW(:), IOW(:), NOCSF(:), IOCSF(:), SCR(:)
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
      IAC=MIN(A0,C0)
      NVERT0=((A0+1)*(C0+1)*(2*B0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
      NDRT0=5*NVERT0
      NDOWN0=4*NVERT0
      NTMP=((NLEV+1)*(NLEV+2))/2
      Call mma_allocate(DRT0,NDRT0,Label='DRT0')
      Call mma_allocate(DOWN0,NDOWN0,Label='DOWN0')
      Call mma_allocate(TMP,NTMP,Label='TMP')
      Call DRT0_MCLR  ! Set up the guga table
     &     (A0,B0,C0,NVERT0,DRT0,DOWN0,NTMP,TMP)

      If ( iPrint.ge.5 )Call PRDRT_MCLR(NVERT0,DRT0,DOWN0)
      Call mma_deallocate(TMP)
*
      LV1RAS=ntRas1
      LV3RAS=LV1RAS+ntRas2
      LM1RAS=2*LV1RAS-nHole1
      LM3RAS=nActEl-nElec3
      Call mma_allocate(V11,NVERT0,Label='V11')
      Call RESTR_MCLR   ! PUT THE RAS CONSTRAINT TO THE DRT TABLE
     &     (NVERT0,DRT0,DOWN0,V11,LV1RAS,LV3RAS,LM1RAS,LM3RAS,NVERT)
*

      NDRT=5*NVERT
      NDOWN=4*NVERT
      Call mma_allocate(DRT,NDRT,Label='DRT')
      Call mma_allocate(DOWN,NDOWN,Label='DOWN')
      Call DRT_MCLR  ! Set up the DRT table used in calculation
     &     (NVERT0,NVERT,DRT0,DOWN0,V11,DRT,DOWN)
      Call mma_deallocate(V11)
      Call mma_deallocate(DOWN0)
      Call mma_deallocate(DRT0)
*
      NDAW=5*NVERT
      Call mma_allocate(DAW,NDAW,Label='DAW')
      Call MKDAW_MCLR(NVERT,DOWN,DAW,iPrint)
*
      NUP=4*NVERT
      NRAW=5*NVERT
      Call mma_allocate(UP,NUP,Label='UP')
      Call mma_allocate(RAW,NRAW,Label='RAW')
      Call MKRAW_MCLR
     &     (NVERT,DOWN,DAW,UP,RAW,iPrint)
*
      NLTV=NLEV+2
      Call mma_allocate(LTV,NLTV,Label='LTV')
      Call MKMID_MCLR
     &     (NVERT,NLEV,DRT,DOWN,DAW,UP,RAW,LTV,
     &      MIDLEV,NMIDV,MIDV1,MIDV2,MXUP,MXDWN,iPrint)
      Call mma_deallocate(LTV)
*
      NIPWLK=1+(MIDLEV-1)/15
      NIPWLK=MAX(NIPWLK,1+(NLEV-MIDLEV-1)/15)
      NNOW=2*NMIDV*nSym
      NIOW=NNOW
      NNOCSF=NMIDV*(nSym**2)
      NIOCSF=NNOCSF
      NSCR=3*(NLEV+1)
      Call mma_allocate(NOW,NNOW,Label='NOW')
      Call mma_allocate(IOW,NIOW,Label='IOW')
      Call mma_allocate(NOCSF,NNOCSF,Label='NOCSF')
      Call mma_allocate(IOCSF,NIOCSF,Label='IOCSF')
      Call mma_allocate(SCR,NSCR,Label='SCR')
      Call MKCOT_MCLR
     &     (nSym,NLEV,NVERT,MIDLEV,NMIDV,MIDV1,MIDV2,NWALK,NIPWLK,
     &      OrbSym,DOWN,NOW,IOW,NCSF,IOCSF,NOCSF,SCR,iPrint)
*
      If ( nConf.ne.NCSF(state_sym).and.(nConf.ne.1) ) then
         Write (6,*) "Set nConf=NCSF(state_sym)"
         Write (6,*)
         nConf=NCSF(state_sym)
      End If
*
      Call mma_deallocate(SCR)
      Call mma_deallocate(IOCSF)
      Call mma_deallocate(NOCSF)
      Call mma_deallocate(IOW)
      Call mma_deallocate(NOW)
      Call mma_deallocate(RAW)
      Call mma_deallocate(UP)
      Call mma_deallocate(DAW)
      Call mma_deallocate(DOWN)
      Call mma_deallocate(DRT)

*
*

      Return
      End
