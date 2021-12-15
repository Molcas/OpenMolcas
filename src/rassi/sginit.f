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
      Subroutine SGInit(nSym,nActEl,iSpin,nRasPrt,nRas,nRasEl,
     &                  iSGStruct)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rassi.fh"
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SGINIT')
      dimension nRas(8,nRasPrt),nRasEl(nRasPrt)
#include "Struct.fh"
      Dimension iSGStruct(nSGSize)
#include "WrkSpc.fh"




      NLEV=NASHT
C Allocate Level to Symmetry table ISm:
      Call GetMem('ISm','Allo','Integer',lISm,nLev)
      ITABS=0
      IT=0
      DO ISYM=1,NSYM
        DO IT=1,NASH(ISYM)
          ITABS=ITABS+1
          ILEV=LEVEL(ITABS)
          IWORK(LISM-1+ILEV)=ISYM
        END DO
      END DO

C Compute size of unrestricted DRT table:
      ib0=ispin-1
      ia0=(nActEl-ib0)/2
      ic0=nLev-ia0-ib0
CTEST      write(*,*)' ia0, ib0, ic0:'
CTEST      write(*,'(8i5)')ia0,ib0,ic0

      iErr=0
      If ((2*ia0+ib0).ne.nActEl) Then
        iErr=1
      Else If((ia0.lt.0).or.(ib0.lt.0).or.(ic0.lt.0)) Then
        iErr=1
      End If
      If(iErr.ne.0) then
        Write(6,*)' RASSI/SGINIT: Impossible input variables.'
        Write(6,*)'   nLev:',nLev
        Write(6,*)' nActEl:',nActEl
        Write(6,*)'  iSpin:',iSpin
        Write(6,*)'Program stops, sorry.'
        CALL ABEND()
      End If

      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
CTEST      write(*,*)' NVERT0:',NVERT0
      NDRT0=5*NVERT0
      NDOWN0=4*NVERT0

C Compute unrestricted DRT tables:
      CALL GETMEM('DRT0  ','ALLO','INTEGER',LDRT0,NDRT0)
      CALL GETMEM('DOWN0 ','ALLO','INTEGER',LDOWN0,NDOWN0)
      NTMP=((NLEV+1)*(NLEV+2))/2
      CALL GETMEM('TMP   ','ALLO','INTEGER',LTMP,NTMP)
      CALL DRT0_RASSI (IA0,IB0,IC0,NVERT0,IWORK(LDRT0),IWORK(LDOWN0),
     &           NTMP,IWORK(LTMP))
CTEST      write(*,*)' SGINIT: Back from DRT0'

C Construct a restricted graph.
      Call GetMem('Lim  ','Allo','Inte',lLim,nLev)
      Do Lev=1,nLev
        IWork(lLim-1+Lev)=0
      End Do
CTEST      write(*,*)' Initialize Work(LLIM):'
CTEST      write(*,'(1x,10i5)')(IWork(llim-1+i),i=1,nlev)
C Fill in the occupation limit table:
      Lev=0
      Do iRO=1,nRasPrt
        Do iSy=1,nSym
          Lev=Lev+nRas(iSy,iRO)
        End Do
        if(Lev.gt.0) IWork(lLim-1+Lev)=nRasEl(iRO)
      End Do
CTEST      write(*,*)' After filling in first values, Work(LLIM):'
CTEST      write(*,'(1x,10i5)')(IWork(llim-1+i),i=1,nlev)
      Call GetMem('NwVer ','Allo','Inte',lNWV,nVert0)
      nVert=nVert0
      Call RmVert(nLev,nVert,IWork(lDRT0),IWork(lDown0),
     &              IWork(lLim),IWork(lNWV))
CTEST      write(*,*)' Back from RMVERT'
      Call GetMem('Lim  ','Free','Inte',lLim,nLev)
      Call GetMem('DRT','Allo','Inte',lDRT,5*nVert)
      Call GetMem('Down','Allo','Inte',lDown,4*nVert)
      Call DRT_RASSI(nVert0,IWork(lDRT0),IWork(lDown0),IWork(lNWV),
     &         nVert,IWork(lDRT),IWork(lDown))
CTEST      write(*,*)' Back from DRT. NVERT=',NVERT
      Call GetMem('NwVer ','Free','Inte',lNWV,NVERT0)
      CALL GETMEM('      ','FREE','Inte',LDRT0,NDRT0)
      CALL GETMEM('      ','FREE','Inte',LDOWN0,NDOWN0)

C Direct Arc Weights table and Level-To-Vertex table:
      Call GetMem('DAW','Allo','Inte',lDAW,5*nVert)
      Call GetMem('LTV','Allo','Inte',lLTV,nLev+2)
      Call MkDAW_RASSI(nLev,nVert,IWork(lDRT),IWork(lDown),IWork(lDAW),
     &           IWork(lLTV))

C Upchain Index table:
      Call GetMem('UP','Allo','Inte',lUp,4*nVert)
C Reverse Arc Weights table:
      Call GetMem('RAW','Allo','Inte',lRAW,5*nVert)
C Modified Arc Weights table:
      Call GetMem('MAW','Allo','Inte',lMAW,4*nVert)
      Call MkMAW(nLev,nVert,IWork(lDown),IWork(lDAW),iWork(lUp),
     &           IWork(lRAW),IWork(lMAW),IWork(lLTV),MidLev)
      MVSta=IWork(lLTV+1+MidLev)
      MVEnd=IWork(lLTV+MidLev)-1
C The DAW, RAW tables are no longer needed:
      Call GetMem('RAW','Free','Inte',lRAW,5*nVert)
      Call GetMem('DAW','Free','Inte',lDAW,5*nVert)
      CALL GETMEM('TMP   ','FREE','INTEGER',LTMP,NTMP)

C Put sizes and addresses in structure iSGStruct:
      iSGStruct(1) =nSym
      iSGStruct(2) =nLev
      iSGStruct(3) =lISm
      iSGStruct(4) =nVert
      iSGStruct(5) =lDRT
      iSGStruct(6) =lDown
      iSGStruct(7) =lUp
      iSGStruct(8) =MidLev
      iSGStruct(9) =MVSta
      iSGStruct(10) =MVEnd
      iSGStruct(11)=lMAW
      iSGStruct(12)=lLTV

      return
      end
