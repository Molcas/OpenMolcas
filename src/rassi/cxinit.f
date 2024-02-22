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
      Subroutine CXInit(SGS,CIS,iXStruct)
      use stdalloc, only: mma_allocate
      use Struct, only:  nXSize, SGStruct, CIStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Dimension iXStruct (nXSize)
#include "WrkSpc.fh"

      nSym   =SGS%nSym
      nLev   =SGS%nLev
      nVert  =SGS%nVert
      MidLev =SGS%MidLev
      MVSta  =SGS%MVSta
      MVEnd  =SGS%MVEnd

C Calculate segment values, and MVL and MVR tables:
      nMidV=1+MVEnd-MVSta
C nIpWlk: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      nIpWlk=1+(MidLev-1)/15
      nIpWlk=MAX(nIpWlk,1+(nLev-MidLev-1)/15)
      Call GetMem('IVR','Allo','Inte',lIVR,2*nVert)
      Call GetMem('MVR','Allo','Inte',lMVR,2*nMidV)
      Call GetMem('MVL','Allo','Inte',lMVL,2*nMidV)
      nSgmnt=26*nVert
      Call GetMem('ISGM','Allo','Inte',lISgm,nSgmnt)
      Call GetMem('VSGM','Allo','Real',lVSgm,nSgmnt)
CTEST      write(*,*)' Calling MKSEG.'
      Call MkSeg(SGS,nLev,nVert,nMidv,
     &        SGS%DRT,SGS%Down,SGS%LTV,
     &        IWork(lIVR),IWork(lMVL),IWork(lMVR),
     &        IWork(lISgm),Work(lVSgm))
CTEST      write(*,*)' Back from MKSEG.'
      CIS%nMidV   =nMidV
      CIS%nIpWlk  = nIpWlk
      iXStruct(8 )=lMVL
      iXStruct(9 )=lMVR

C Various offset tables:
      nNOW=2*nMidV*nSym
CUNUSED      nIOW=nNOW
      Call mma_allocate(CIS%NOW,nNOW,Label='CIS%NOW')
      Call mma_allocate(CIS%IOW,nNOW,Label='CIS%IOW')
      MxEO=(nLev*(nLev+5))/2
      nNOCP=MxEO*nMidV*nSym
      nIOCP=nNOCP
      nNRL=(1+MxEO)*nVert*nSym
      Call GetMem('NOCP','Allo','Inte',lNOCP,nNOCP)
      Call GetMem('IOCP','Allo','Inte',lIOCP,nIOCP)
      Call mma_allocate(CIS%NCSF,nSym,Label='CIS%NCSF')
      Call GetMem('NRL','Allo','Inte',lNRL,nNRL)
      nNOCSF=nMidV*(nSym**2)
      nIOCSF=nNOCSF

      Call GetMem('NOCSF','Allo','Inte',lNOCSF,nNOCSF)
      Call GetMem('IOCSF','Allo','Inte',lIOCSF,nIOCSF)
      CIS%lNOCSF  =lNOCSF
      CIS%lIOCSF  =lIOCSF
      iXStruct(1)=MxEO
      iXStruct(2)=lNOCP
      iXStruct(3)=lIOCP
      Call NrCoup(SGS,CIS,iXStruct,
     &         nVert,nMidV,MxEO,SGS%ISm,SGS%DRT,
     &         IWork(lISgm),CIS%NOW,CIS%IOW,IWork(lNOCP),
     &         IWork(lIOCP),IWork(lNOCSF),IWork(lIOCSF),
     &         CIS%NCSF,IWork(lNRL),IWork(lMVL),IWork(lMVR))
CTEST      write(*,*)' Back from NRCOUP.'
      Call GetMem('NRL','Free','Inte',lNRL,nNRL)
C Computed in NrCoup:
      nWalk=CIS%nWalk
      nICoup=IXSTRUCT(4)

      nICase=nWalk*nIpWlk
CTEST      write(*,*)' NWALK:',NWALK
      Call GetMem('ICASE','Allo','Inte',lICase,nICase)
      nnICoup=3*nICoup
      Call GetMem('ICoup','Allo','Inte',lICoup,nnICoup)
      nVMax=5000
      Call GetMem('VTabTmp','Allo','Real',lVTabTmp,nVMax)
      nILNDW=nWalk
      Call GetMem('iLndw','Allo','Inte',liLndw,niLndw)
      nScr=7*(nLev+1)
      Call GetMem('SCR','Allo','Inte',lScr,nScr)
      Call GetMem('VAL','Allo','Real',lVal,nLev+1)
      CIS%lICase=lICase
      iXStruct(5)=lICoup
      iXStruct(6)=nVMax
      iXStruct(7)=lVTabTmp
      nVTab=nVMax
      lVTab=lVTabTmp
      Call MkCoup(nLev,SGS%Ism,nVert,MidLev,nMidV,MVSta,MVEnd,
     &            MxEO,nICoup,nWalk,nICase,nVTab,
     &            IWork(lIVR),SGS%MAW,IWork(lISGM),
     &            WORK(lVSGM),CIS%NOW,CIS%IOW,IWork(lNOCP),
     &     IWork(lIOCP),IWork(lILNDW),IWork(lICASE),IWork(lICOUP),
     &     WORK(lVTAB),IWork(lSCR),WORK(lVAL))
C iXStruct(10)..iXStruct(14) are set in MkCoup

C nVTab has now been updated to the true size. Allocate final array:
      Call GetMem('VTab','Allo','Real',lVtab,nVTab)
      iXStruct(6)=nVTab
      iXStruct(7)=lVTab
      call dcopy_(nVTab,Work(lVTabTmp),1,Work(lVTab),1)
      Call GetMem('VTabTmp','Free','Real',lVTabTmp,nVMax)
      Call GetMem('iLndw','Free','Inte',liLndw,niLndw)
      Call GetMem('SCR','Free','Inte',lScr,nScr)
      Call GetMem('VAL','Free','Real',lVal,nLev+1)
      Call GetMem('ISGM','Free','Inte',lISgm,nSgmnt)
      Call GetMem('VSGM','Free','Real',lVSgm,nSgmnt)
      Call GetMem('IVR','Free','Inte',lIVR,2*nVert)

      end
