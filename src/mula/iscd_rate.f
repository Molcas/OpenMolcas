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
* Copyright (C) 2009, Giovanni Ghigo                                   *
************************************************************************
C!-----------------------------------------------------------------------!
C!
C!  InterSystem Crossing rate
C!  Author: Giovanni Ghigo
C!          Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
C!          June 2009
C!-----------------------------------------------------------------------!
      Subroutine ISCD_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,nYes,
     &           dMinWind,lBatch,nBatch,leftBatch,nIndex,VibWind2,
     &           lNMAT0,lNMAT,lNINC,lNDEC,lnTabDim,nnTabDim,
     &           C1,C2,W1,W2,det0,det1,det2,
     &           C,W,r01,r02,r00,
     &           m_max,n_max,max_dip,nnsiz,FC00,FCWind2,dRho,
     &           mTabDim,mMat,mInc,mDec,
     &                   nMat,nInc,nDec)
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
#include "Constants_mula.fh"
#include "inout.fh"
#include "WrkSpc.fh"
#include "io_mula.fh"
      Integer nIndex(3,0:maxMax_n)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nosc),W1(nOsc,nOsc),W2(nOsc,nOsc),
     &       C(nOsc,nOsc),W(nOsc,nOsc)
      Real*8 r01(nOsc),r02(nOsc),r00(nOsc), det0,det1,det2, FC00
      Real*8 FCWind2(nYes)
      Integer VibWind2(nYes), nnTabDim(0:lnTabDim)
      Integer mMat(0:mTabDim,nOsc) , nMat(nOsc,lBatch)
      Integer mInc(0:mTabDim,nOsc) , nInc(nOsc,lBatch)
      Integer mDec(0:mTabDim,nOsc) , nDec(nOsc,lBatch)

      Call TabDim2_drv(m_max,nosc,nvTabDim)
      Call TabDim2_drv(n_max,nosc,nvTabDim)
      max_nOrd = nvTabDim-1
      Call TabDim2_drv(m_max,nosc,nvTabDim)
      m_max_ord  = nvTabDim-1
      Call TabDim2_drv(min(n_max,m_max+1),nosc,nvTabDim)
      mx_max_ord = nvTabDim-1
      Call TabDim2_drv(min(m_max,n_max+1),nosc,nvTabDim)
      nx_max_ord = nvTabDim-1
      Call TabDim2_drv(m_max-1,nosc,nvTabDim)
      max_mInc   = nvTabDim-1
      Call TabDim2_drv(n_max-1,nosc,nvTabDim)
      max_nInc   = nvTabDim-1
      Call TabDim2_drv(n_max,nosc,nvTabDim)
      n_max_ord  = nvTabDim-1

      mx_max_ord = 0 ! CGGn
      If (iPrint.GE.3) Write(6,*) ' Memory allocated for U matrix:',
     &           (n_max_ord+1)*(mx_max_ord+1),          ' words,  ',
     &         8*(n_max_ord+1)*(mx_max_ord+1)/1024/1024,' MB.     '
      Call XFlush(6)
      Call GetMem('L','Allo','Real', ipL,(m_max_ord+1)*(nx_max_ord+1))
      Call GetMem('U','Allo','Real',ipU,(n_max_ord+1)*(mx_max_ord+1))
      Call GetMem('alpha1','Allo','Real',ipAlpha1,nOsc*nOsc)
      Call GetMem('alpha2','Allo','Real',ipAlpha2,nOsc*nOsc)
      Call GetMem('beta','Allo','Real',ipBeta,nOsc*nOsc)
      Call GetMem('MAT0','Allo','Inte',ipnMat0,nOsc)

CGGt -------------------------------------------------------------------
c      Write(6,*)'     lnTabDim+1=',lnTabDim+1,':'                 ! CGGt
c      Do i = 0, lnTabDim                                          ! CGGt
c      iIndex0 = nnTabDim(i)                                       ! CGGt
c      Call iDaFile(lNMAT0,2,iWork(ipnMat0),nOsc,iIndex0)          ! CGGt
c      Write(6,*) i,' read at',nnTabDim(i),'  M:',                 ! CGGt
c     &                              (iWork(ipnMat0+j),j=0,nOsc-1) ! CGGt
c      EndDo                                                       ! CGGt
c      Write(6,*)'-----------------------------------------------' ! CGGt
CGGt -------------------------------------------------------------------
c      Call GetMem('Test_2','LIST','INTE',iDum,iDum)               ! CGGt
c      Call XFlush(6)                                              ! CGGt
      Call ISCD_FCval(iPrint,   iMaxYes,  lnTabDim,      nnTabDim,
     &  lNMAT0,    lNMAT,      lNINC,     lNDEC,
     &  lBatch,    nBatch,     leftBatch, nIndex,
     &  C1,        W1,         det1,      r01,           C2,
     &  W2,        det2,       r02,       m_max_ord,
     &  n_max_ord, mx_max_ord, max_mInc,  max_nInc,      nx_max_ord,
     &  mMat,      nMat,       mInc,      nInc,          mDec,
     &  nDec,      C,          W,         det0,          r00,
     &  Work(ipL), Work(ipU),  FC00,      Work(ipAlpha1),Work(ipAlpha2),
     &  Work(ipBeta), nOsc, nnsiz, iMx_nOrd,
     &  nYes,      VibWind2,   FCWind2,   iWork(ipnMat0) )

      const = 2.0d0*rpi/5.309d-12
      If (iPrint.GE.4) then
        Write(6,*)
        Write(6,*)'  const =',const
        Write(6,*)'  dRho/cm =', dRho/HarToRcm
        Write(6,*)'  const*dRho=',const*dRho/HarToRcm
      EndIf
      const = const*dRho/HarToRcm

      dSum = 0.0d0
      Do ii = 1, nYes
        dSum = dSum + FCWind2(ii)**2
      EndDo

      dSoc = 0.0d0
      Do ii = 1, 3
        dSOC = dSOC + TranDip(ii)**2
      EndDo

      dRate = const * dSum * dSoc / dMinWind
      dLT = 1.0d0/dRate

      If (iPrint.GE.3) then
        Write(6,*) '  Sum of squares of FC factors =',dSum
        Write(6,*) '  Root-square of the sum =', SQRT(dSum)
        Write(6,*) '  dSOC =', dSOC
      EndIf

      If (iPrint.GE.1) then
        Write(6,*)
        Write(6,*) ' InterSystem Crossing rate constant: '
        Write(6,*) ' ===================================='
        Write(6,'(a,e10.2,a)') '  ISC Rate Constant  ',dRate,' sec-1'
        Write(6,'(a,e10.2,a)') '  Lifetime           ',dLT,' sec'
        dLT = dLT * 1.0d3
        If (dLT.GT.1.0d0 .and. dLT.LE.1.0d3)
     &  Write(6,'(a19,f5.1,a)') ' ',dLT,' msec'
        dLT = dLT * 1.0d3
        If (dLT.GT.1.0d0 .and. dLT.LE.1.0d3)
     &  Write(6,'(a19,f5.1,a)') ' ',dLT,' microsec'
        dLT = dLT * 1.0d3
        If (dLT.GT.1.0d0 .and. dLT.LE.1.0d3)
     &  Write(6,'(a19,f5.1,a)') ' ',dLT,' nsec'
        dLT = dLT * 1.0d3
        If (dLT.GT.1.0d0 .and. dLT.LE.1.0d3)
     &  Write(6,'(a19,f5.1,a)') ' ',dLT,' psec'
        dLT = dLT * 1.0d3
        If (dLT.GT.1.0d0 .and. dLT.LE.1.0d3)
     &  Write(6,'(a19,f5.1,a)') ' ',dLT,' fsec'
        Write(6,*) ' ------------------------------------'
        Write(6,*)
        Write(6,*)
        Call XFlush(6)
      EndIf

      Call GetMem('MAT0','Free','Inte',ipnMat0,nOsc)
      Call GetMem('beta','Free','Real',ipBeta,nOsc*nOsc)
      Call GetMem('alpha2','Free','Real',ipAlpha2,nOsc*nOsc)
      Call GetMem('alpha1','Free','Real',ipAlpha1,nOsc*nOsc)
      Call GetMem('U','Free','Free',ipU,(n_max_ord+1)*(mx_max_ord+1))
      Call GetMem('L','Free','Real', ipL,(m_max_ord+1)*(nx_max_ord+1))

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(max_dip)
      End
