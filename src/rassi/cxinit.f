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
      Subroutine CXInit(SGS,CIS,EXS)
      use stdalloc, only: mma_allocate
      use Struct, only:  SGStruct, CIStruct, EXStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS
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
      Call mma_allocate(EXS%MVR,2*nMidV)
      Call mma_allocate(EXS%MVL,2*nMidV)
      nSgmnt=26*nVert
      Call GetMem('ISGM','Allo','Inte',lISgm,nSgmnt)
      Call GetMem('VSGM','Allo','Real',lVSgm,nSgmnt)
      Call MkSeg(SGS,nLev,nVert,nMidv,
     &        SGS%DRT,SGS%Down,SGS%LTV,
     &        IWork(lIVR),EXS%MVL,EXS%MVR,
     &        IWork(lISgm),Work(lVSgm))
      CIS%nMidV   =nMidV
      CIS%nIpWlk  = nIpWlk

C Various offset tables:
      nNOW=2*nMidV*nSym
      Call mma_allocate(CIS%NOW,nNOW,Label='CIS%NOW')
      Call mma_allocate(CIS%IOW,nNOW,Label='CIS%IOW')
      MxEO=(nLev*(nLev+5))/2
      nNOCP=MxEO*nMidV*nSym
      nIOCP=nNOCP
      nNRL=(1+MxEO)*nVert*nSym
      Call mma_allocate(EXS%NOCP,nNOCP,Label='EXS%NOCP')
      Call mma_allocate(EXS%IOCP,nIOCP,Label='EXS%IOCP')
      Call mma_allocate(CIS%NCSF,nSym,Label='CIS%NCSF')
      Call GetMem('NRL','Allo','Inte',lNRL,nNRL)
      nNOCSF=nMidV*(nSym**2)
      nIOCSF=nNOCSF

      Call mma_allocate(CIS%NOCSF,nNOCSF,Label='CIS%NOCSF')
      Call mma_allocate(CIS%IOCSF,nIOCSF,Label='CIS%IOCSF')
      EXS%MxEO =MxEO
      Call NrCoup(SGS,CIS,EXS,
     &         nVert,nMidV,MxEO,SGS%ISm,SGS%DRT,
     &         IWork(lISgm),CIS%NOW,CIS%IOW,EXS%NOCP,
     &         EXS%IOCP,CIS%NOCSF,CIS%IOCSF,
     &         CIS%NCSF,IWork(lNRL),EXS%MVL,EXS%MVR)
      Call GetMem('NRL','Free','Inte',lNRL,nNRL)
C Computed in NrCoup:
      nWalk=CIS%nWalk
      nICoup=EXS%nICoup

      nICase=nWalk*nIpWlk
      Call mma_allocate(CIS%ICase,nICase,Label='CIS%ICase')
      nnICoup=3*nICoup
      Call mma_allocate(EXS%ICoup,nnICoup,Label='EXS%ICoup')
      nVMax=5000
      Call GetMem('VTabTmp','Allo','Real',lVTabTmp,nVMax)
      nILNDW=nWalk
      Call GetMem('iLndw','Allo','Inte',liLndw,niLndw)
      nScr=7*(nLev+1)
      Call GetMem('SCR','Allo','Inte',lScr,nScr)
      Call GetMem('VAL','Allo','Real',lVal,nLev+1)
      EXS%nVTab =nVMax
      nVTab=nVMax
      lVTab=lVTabTmp
      Call MkCoup(nLev,SGS%Ism,nVert,MidLev,nMidV,MVSta,MVEnd,
     &            MxEO,nICoup,nWalk,nICase,nVTab,
     &            IWork(lIVR),SGS%MAW,IWork(lISGM),
     &            WORK(lVSGM),CIS%NOW,CIS%IOW,EXS%NOCP,
     &            EXS%IOCP,IWork(lILNDW),CIS%ICase,EXS%ICOUP,
     &            WORK(lVTAB),IWork(lSCR),WORK(lVAL))

C nVTab has now been updated to the true size. Allocate final array:
      Call mma_allocate(EXS%Vtab,nVTab,Label='EXS%VTab')
      EXS%nVTab=nVTab
      call dcopy_(nVTab,Work(lVTabTmp),1,EXS%VTab,1)
      Call GetMem('VTabTmp','Free','Real',lVTabTmp,nVMax)
      Call GetMem('iLndw','Free','Inte',liLndw,niLndw)
      Call GetMem('SCR','Free','Inte',lScr,nScr)
      Call GetMem('VAL','Free','Real',lVal,nLev+1)
      Call GetMem('ISGM','Free','Inte',lISgm,nSgmnt)
      Call GetMem('VSGM','Free','Real',lVSgm,nSgmnt)
      Call GetMem('IVR','Free','Inte',lIVR,2*nVert)

      end Subroutine CXInit
