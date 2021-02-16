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
C!  InterSystem Crossing rate evaluation: "Engine" routines
C!  Author: Giovanni Ghigo
C!          Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
C!          28 Dec-08 - 06 Jan-09
C!
      Subroutine ISC_Rho(iPrint,nOsc,new_n_max,dRho,energy1,  energy2,
     &                     minQ, dMinWind, nMaxQ, harmfreq1,harmfreq2)
C!
C!      Calculate State Density  dRho  GG 30-Dec-08
C!      Formula (86) taken from  M. Bixon, J. Jortner  JCP,48,715 (1969)
C!
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
#include "Constants_mula.fh"
      Real*8 energy1, energy2, harmfreq1(nOsc), harmfreq2(nOsc)
      Real*8 GE1, GE2, dMinWind, dMinWind0
      Integer nMaxQ(nOsc)
C!
      If (iPrint.GE.2) then
        Write(6,*)
        Write(6,*) ' State Density data:                         '
        Write(6,*) ' ============================================'
      EndIf
C!
      dMinWind0 = dMinWind
      If (dMinWind.EQ.0.0d0) dMinWind = 1.0d0
      minQ = 0
      dRho = 2.0d0 / rpi
      dRho = SQRT(dRho / 1.0d0*nOsc)
      dRho = dRho * (1.0d0 - 1.0d0/(12.0d0*nOsc))
      avFreq = 0.0d0
      avFreqSq = 0.0d0
      dLambda = 1.0d0
      dZPE1 = 0.0d0
      dZPE2 = 0.0d0
      dMinFreq2 = 1.0d99
      dMaxFreq2 = 0.0d0
      Do iOsc = 1,nOsc
        avFreq   = avFreq + harmfreq2(iOsc)
        avFreqSq = avFreqSq + (harmfreq2(iOsc))**2
        dZPE1 = dZPE1 + 0.5d0*harmfreq1(iOsc)
        dZPE2 = dZPE2 + 0.5d0*harmfreq2(iOsc)
        dMinFreq2 = min( dMinFreq2 , harmfreq2(iOsc) )
        dMaxFreq2 = Max( dMaxFreq2 , harmfreq2(iOsc) )
      EndDo
      GE1 = energy1 + dZPE1
      GE2 = energy2 + dZPE2
      T0 = ABS(GE2 - GE1)
      avFreq = avFreq / nOsc
      avFreqSq = avFreqSq / nOsc

      Do jOsc = 1,nOsc
        dLambda = dLambda *  harmfreq2(jOsc)/avFreq
      EndDo
      new_n_max = INT(0.5d0 + (GE1-GE2) / dMinFreq2)
      dAlpha = avFreqSq/(avFreq**2)
      dBeta = ( (nOsc-1) * (nOsc-2) * dAlpha -nOsc**2 ) / (6*nOsc)
      dEtha = ABS(GE1-GE2) / dZPE2
      dRho = dRho / avFreq
      dRho = dRho * dAlpha
      dRho = dRho / (1+dEtha)
      dDn = ( 1.0d0 + (2.0d0/dEtha))
      dDn = dDn ** (dEtha/2.0d0)
      dDn = dDn * ( 1.0d0 + (dEtha/2.d0) )
      dDn = dDn ** nOsc
      dRho = dRho * dDn
      dFE = ( 1.0d0 + dEtha )**2
      dFE = ( 1.0d0 - 1.0d0/dFE )**dBeta
      dRho = dRho * dFE

      Do jOsc = 1,nOsc
        nMaxQ(jOsc) = INT(0.5d0+(T0+dMinWind/dRho) / harmfreq2(jOsc) )
      EndDo
      minQ = INT(0.5d0+(T0-dMinWind/dRho) / dMaxFreq2 )
      If (minQ.LT.0) then
        Write(6,*)
        Write(6,*) ' ***** ERROR ******'
        write(6,*) ' Window too large !'
        Write(6,*) ' ******************'
        Call Quit_OnUserError()
      EndIf

      If (iPrint.GE.2) then
        Write(6,'(a,f11.6,a)')  '  T_0  = ',T0,         ' (au)  '
        Write(6,'(a,f11.3,a)')  '  T_0  = ',T0*HarToRcm,' (cm-1)'
        Write(6,'(a,f11.3,a)')  '  T_0  = ',T0*27.2114, ' (eV)  '
        Write(6,'(a,d14.3,a)')  '  State Density (dRho) = ',
     &                 dRho,' (au-1)'
        Write(6,'(a,g14.3,a)')  '  State Density (dRho) = ',
     &                 dRho/HarToRcm,' (cm)'
        Write(6,'(a,g17.9,a)') '  1/dRho = ',
     &                       HarToRcm/dRho,' (cm-1)'
        Write(6,'(a,f7.3,a)')  '  Expansion factor =',dMinWind
        Write(6,'(a,g17.9,a)') '  Window = (+/-)',
     &        0.5d0*dMinWind*HarToRcm/dRho,' (cm-1)'
      EndIf
      If (iPrint.GE.3) then
        Write(6,*) ' Maximum quantum numbers:',(nMaxQ(i),i=1,nOsc)
        Write(6,*) ' Minimum quantum number: ',minQ
        Write(6,*) ' Suggested n_max (new_n_max)=',new_n_max
        Write(6,*)
      EndIf
      Call XFlush(6)
      dMinWind = dMinWind0
C!
      Return
      End


      Subroutine ISC_Ene(iPrint,nOsc,max_nOrd,nYes,nMat,nTabDim,
     & GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,
     & dMinWind,dRho, EneMat, lVec, lTVec)
cC!
cC!    Calculate Energy of Levels  GG 30-Dec-08 - 08-Jan-09
cC!
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
#include "Constants_mula.fh"
#include "WrkSpc.fh"
      Real*8 GE1, GE2, harmfreq1(nOsc), harmfreq2(nOsc)
      Real*8 x_anharm1(nOsc,nOsc),x_anharm2(nOsc,nOsc)
      Real*8 dMinWind, dRho, dWlow, dWup
      Real*8 dEne, EneMat(0:max_nOrd)
      Integer nMat(0:nTabDim,nOsc), lVec(0:nTabDim), lTVec(0:nTabDim)
      Logical lUpdate

      If (dMinWind.EQ.0.0d0) then
        lUpDate=.True.
        dMinWind = 1.0d0
      else
        lUpDate=.False.
      EndIf
      Do iOrd = 0, max_nOrd
        lTVec(iOrd) = lVec(iOrd)
      EndDo
C!
C!    Energy calculation
C!
      If (iPrint.GE.4) then
        Write(6,*)
        Write(6,*) ' States in the preliminar window :'
        If (nOsc.LE.24) then
          Write(6,'(a,108a)') '  ',('=',i=1,108)
          Write(6,*)'     jOrd    ene/au    ene/cm-1 '//
     &    'Vibrational quantum numbers'
          Write(6,'(a,108a)') '  ',('-',i=1,108)
        else
          Write(6,'(a,36a)') '  ',('=',i=1,36)
          Write(6,*)'        #    jOrd   ene/au      ene/cm-1 '
          Write(6,'(a,36a)') '  ',('-',i=1,36)
        EndIf
        Call XFlush(6)
      EndIf
C!
      Call GetMem('level1','Allo','Inte',iplevel1,nOsc)
      Call GetMem('level2','Allo','Inte',iplevel2,nOsc)
      Do iv=1,nOsc
        iWork(iplevel1+iv-1) = 0
      EndDo
      Do iOrd = 0, max_nOrd
        If (lVec(iOrd).EQ.1) then
          Do iv=1,nOsc
            iWork(iplevel2+iv-1) = nMat(iOrd,iv)
          EndDo
          l_harm=nOsc
          Call TransEnergy(
     &    GE1,x_anharm1,harmfreq1,iWork(iplevel1),
     &    GE2,x_anharm2,harmfreq2,iWork(iplevel2),
     &    dEne,l_harm)
          EneMat(iOrd) = dEne
          If (iPrint.GE.4) then
            If (nOsc.LE.24) then
              loc_n_max = 0
              Do j=1,nOsc
                loc_n_max = loc_n_max + nMat(iOrd,j)
              EndDo
              Write(6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,
     &            dEne*HarToRcm,loc_n_max,': ',(nMat(iOrd,j),j=1,nOsc)
            else
               Write(6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,
     &            dEne*HarToRcm,loc_n_max
            EndIf
          EndIf
        EndIf
      End Do
C!
C!    Energy selection
C!
      If (iPrint.GE.3) then
        Write(6,*)
        Write(6,*) ' States in the window :'
        If (nOsc.LE.24) then
          Write(6,'(a,108a)') '  ',('=',i=1,108)
          Write(6,*)'     jOrd    ene/au    ene/cm-1 '//
     &    'Vibrational quantum numbers'
          Write(6,'(a,108a)') '  ',('-',i=1,108)
        else
          Write(6,'(a,36a)') '  ',('=',i=1,36)
          Write(6,*)'        #    jOrd   ene/au      ene/cm-1 '
          Write(6,'(a,36a)') '  ',('-',i=1,36)
        EndIf
        Call XFlush(6)
      EndIf
C!
      nYes_start = nYes
 100  Continue
      dWlow = 0.5d0*dMinWind/dRho
      dWup  = dWlow
      Do iOrd = 0,max_nOrd
        lVec(iOrd) = lTVec(iOrd)
        If (lVec(iOrd).EQ.1) then
          dEne = EneMat(iOrd)
          If (dEne.LT.-dWlow .or. dEne.GT.dWup) then
            lVec(iOrd) = 0
            nYes = nYes - 1
          else
            lVec(iOrd) = 1
            If (iPrint.GE.3) then
              If (nOsc.LE.24) then
                loc_n_max = 0
                Do j=1,nOsc
                  loc_n_max = loc_n_max + nMat(iOrd,j)
                EndDo
                Write(6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,
     &              dEne*HarToRcm,loc_n_max,': ',(nMat(iOrd,j),j=1,nOsc)
              else
                Write(6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,
     &              dEne*HarToRcm,loc_n_max
              EndIf
            EndIf
          EndIf
        EndIf
      EndDo
      If (nYes.LT.1 .and. lUpDate) then
        dMinWind = dMinWind + 1.0d0
        nYes = nYes_start
        GoTo 100
      EndIf
C!
      Call GetMem('level2','Free','Inte',iplevel2,nOsc)
      Call GetMem('level1','Free','Inte',iplevel1,nOsc)

      If (iPrint.GE.3) then
        If (nOsc.LE.30) Write(6,'(a,108a)') '  ',('-',i=1,108)
        If (nOsc.GT.30) Write(6,'(a,36a)') '  ',('-',i=1,36)
        Write(6,'(a,f12.9,a,f12.9,a)')'  Window: ',
     &   -dWlow,         ' / ',dWup,         ' (au)'
        Write(6,'(a,f12.6,a,f12.6,a)')'  Window: ',
     &   -dWlow*HarToRcm,' / ',dWup*HarToRcm,' (cm-1)'
      EndIf
      If (iPrint.GE.2) then
        Write(6,*) ' Final number of States=',nYes
      EndIf
      If (dMinWind.GT.1.0d0 .and. lUpDate .and. iPrint.GE.1) then
        Write(6,*)
        Write(6,*) ' *** Warning: Expansion factor has been set to ',
     &                                                     dMinWind
        Write(6,*)
      EndIf
      Call XFlush(6)
      Return
      End


      Subroutine ISC_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,
     &                    nYes,dMinWind,     VibWind2,
     &                    C1,C2,W1,W2,det0,det1,det2,
     &                    C,W,r01,r02,r00,
     &                    mTabDim,mMat,nTabDim,nMat,
     &                    mInc,mDec,nInc,nDec,
     &                    m_max,n_max,max_dip,nnsiz,FC00,FCWind2,dRho)
C!
C!    Estimate ISC rate  GG 30-Dec-08
C!
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
#include "Constants_mula.fh"
#include "inout.fh"
#include "WrkSpc.fh"
      Real*8 C1(nOsc,nOsc),C2(nOsc,nosc),W1(nOsc,nOsc),W2(nOsc,nOsc),
     &       C(nOsc,nOsc),W(nOsc,nOsc)
      Real*8 r01(nOsc),r02(nOsc),r00(nOsc), det0,det1,det2, FC00
      Real*8 FCWind2(nYes)
      Integer VibWind2(nYes)
      Integer mMat(0:mTabDim,nOsc), nMat(0:nTabDim,nOsc)
      Integer mInc(0:mTabDim,nOsc), nInc(0:iMaxYes,nOsc)
      Integer mDec(0:mTabDim,nOsc), nDec(0:iMaxYes,nOsc)

      Call TabDim2_drv(m_max,nosc,nvTabDim)
      max_mOrd = nvTabDim-1
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

c      Call GetMem('Test_2','LIST','INTE',iDum,iDum)             ! CGGt
c      Call XFlush(6)                                            ! CGGt
      Call ISC_FCval(iPrint,   iMaxYes,  nTabDim,
     &  C1,        W1,         det1,     r01,      C2,
     &  W2,        det2,       r02,                m_max_ord,
     &  n_max_ord, mx_max_ord, max_mInc, max_nInc, nx_max_ord,
     &  mMat,      nMat,       mInc,     nInc,     mDec,
     &  nDec,      C,          W,        det0,     r00,
     &  Work(ipL), Work(ipU),  FC00,     Work(ipAlpha1),Work(ipAlpha2),
     &  Work(ipBeta), nOsc, nnsiz, iMx_nOrd,
     &  nYes, VibWind2, FCWind2)

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

      Call GetMem('beta','Free','Real',ipBeta,nOsc*nOsc)
      Call GetMem('alpha2','Free','Real',ipAlpha2,nOsc*nOsc)
      Call GetMem('alpha1','Free','Real',ipAlpha1,nOsc*nOsc)
      Call GetMem('U','Free','Free',ipU,(n_max_ord+1)*(mx_max_ord+1))
      Call GetMem('L','Free','Real', ipL,(m_max_ord+1)*(nx_max_ord+1))

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(max_dip)
      End
