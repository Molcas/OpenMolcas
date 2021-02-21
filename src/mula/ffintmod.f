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
* Copyright (C) 1998, Anders Bernhardsson                              *
*               1998, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
c       Module ForceFieldIntMod
C!
C!  Purpose:
C!    Calculates transition intensities using the double harmonic
C!    approximation.
C!
C!  Written by:
C!    Anders Bernhardsson & Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1998.
C!
C!-----------------------------------------------------------------------!

Cvv       Private
C!
c       Contains
C!



C!-----------------------------------------------------------------------!
C!
      Subroutine IntForceField(IntensityMat,TermMat,T0,max_term,FC00,
     &   C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,r00,
     &   m_max,n_max,max_dip,
     &   Trandip,TranDipGrad,
     &   harmfreq1,x_anharm1,harmfreq2,x_anharm2,
     &   mMat,mInc,mDec,nMat,nInc,nDec,OscStr,nnsiz,
     &   max_mOrd,max_nOrd,nDimTot,nOsc)
C!
C!  Purpose:
C!    Calculates the intensities of the different transitions between
C!    the two surfaces.
C!
C!  Input:
C!    FC           : Real*8 two dimensional array -
C!                   Franck-Condon factors.
C!    T0           : Real*8 variable - energy difference between
C!                   the two states.
C!    max_term     : Integer - maximum order of the transition dipole terms.
C!    W1,W2        : Real*8 two dimensional arrays - eigenvectors
C!                   scaled by the square root of the eigenvalues. Harmonic
C!                   approximation.
C!    C1,C2        : Real*8 two dimensional arrays - inverses
C!                   of W1 and W2.
C!    det1,det2    : Real*8 variables - determinants of C1 and C2.
C!    r01,r02      : Real*8 arrays - coordinates of the two
C!                   oscillators.
C!    r00          : Real*8 array - coordinates of intermediate
C!                   oscillator.
C!    m_max,n_max  : Integer variables - maximum quanta.
C!    max_dip      : Integer variable - maximum order of transition dipole.
C!    MatEl        : Logical
C!    m_plot,
C!    n_plot       : Integer array - transitions wanted in output.
C!    TermMat      : Real*8 two dimensional array - energies
C!                   of transitions.
C!
C!  Output:
C!    IntensityMat : Real*8 two dimensional array - intensities
C!                   of the transitions.
C!
C!  Uses:
C!    Constants
C!    VibMod
C!    TabMod
C!    FCMod
C!    MatElMod
C!
c       Use VibMod
c       Use TabMod
c       Use FCMod
c       Use MatElMod
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
      Real*8 IntensityMat (0:max_mOrd,0:max_nOrd)
      Real*8 TermMat (0:max_mOrd,0:max_nOrd)
      Real*8 C1(nosc,nosc),C2(nosc,nosc),
     &  W1(nosc,nosc),W2(nosc,nosc),C(nosc,nosc),W (nosc,nosc)
      Real*8 r01(nosc),r02(nosc),r00(nosc)
      Real*8 TranDip (3)
      Real*8 TranDipGrad (3,nosc)
      Real*8 harmfreq1(nosc),harmfreq2  (nosc)
      Real*8 x_anharm1(nosc,nosc),x_anharm2 (nosc,nosc)
      Integer mMat(0:mdim1,mdim2) ,mInc(0:mdim1,mdim2) ,
     &  mDec (0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),
     &  nDec (0:ndim1,ndim2)
      Logical    OscStr
      Integer nvTabDim
#include "WrkSpc.fh"
C!
C!---- Initialize.
      Call TabDim2_drv(m_max,nosc,nvTabDim)
      max_mOrd = nvTabDim-1
      Call TabDim2_drv(n_max,nosc,nvTabDim)
      max_nOrd = nvTabDim-1
      Call GetMem('FC2','Allo','Real',ipFC2,
     &  (max_mOrd+1)*(max_nOrd+1)*4)
C!
      Call SetUpHarmDip(Work(ipFC2),max_dip,m_max,n_max,mMat,mInc,
     &   mDec,nMat,nInc,nDec,
     &   C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,r00,
     &   TranDip,TranDipGrad,FC00,nnsiz,
     &   max_mOrd, max_nOrd, nOsc)
C!
C!---- Calculate intensities with Boltzmann weighting of hotband intensity.
      const1 = (2.0d0/3.0d0)
      Call GetMem('level1','Allo','Inte',iplevel1,nOsc)
      Call GetMem('level2','Allo','Inte',iplevel2,nOsc)

      If ( max_nOrd.gt.max_mOrd ) Then
      l_FreqDiffMat=max_mOrd+1
      Call GetMem('FreqDiffMat','Allo','Real',
     &                ipFreqDiffMat,l_FreqDiffMat)

      do iv=1,nOsc
      iWork(iplevel1+iv-1) = mMat(0,iv)
      enddo
      Do iOrd = 0,max_mOrd
      do iv=1,nOsc
      iWork(iplevel2+iv-1) = mMat(iOrd,iv)
      enddo
      l_harm=nOsc
      Call TransEnergy(0.0d0,x_anharm1,
     &       harmfreq1,iWork(iplevel1),
     &       0.0d0,x_anharm1,harmfreq1,iWork(iplevel2),
     &       Work(ipFreqDiffMat+iOrd),l_harm)
      End Do
      Do jOrd = 0,max_nOrd
      Do iOrd = 0,max_mOrd
      dE = Work(ipFreqDiffMat+iOrd)*HarToaJ*1.0d-18
      const2 = const1*exp(-dE/(1.38066d-23*Temperature))
      IntensityMat(iOrd,jOrd) =
     &       const2*Abs(TermMat(iOrd,jOrd))*
     &          (Work(ipFC2+iOrd+(max_mOrd+1)*
     &          (jOrd +(max_nOrd+1)))**2+
     &           Work(ipFC2+iOrd+(max_mOrd+1)*
     &          (jOrd +(max_nOrd+1)*2))**2+
     &           Work(ipFC2+iOrd+(max_mOrd+1)*
     &          (jOrd +(max_nOrd+1)*3))**2)
      If ( .not.OscStr ) Then
      IntensityMat(iOrd,jOrd) = 32.13002d9*
     &             const2*(TermMat(iOrd,jOrd)**2)*
     &             IntensityMat(iOrd,jOrd)
      End If
      End Do
      End Do
      Else
      l_FreqDiffMat=max_nOrd+1
      Call GetMem('FreqDiffMat','Allo','Real',
     &            ipFreqDiffMat,l_FreqDiffMat)
      do iv=1,nOsc
      iWork(iplevel1+iv-1) = nMat(0,iv)
      enddo
      Do iOrd = 0,max_nOrd
      Do iv=1,nOsc
      iWork(iplevel2+iv-1) = nMat(iOrd,iv)
      enddo
      l_harm=nOsc
      Call TransEnergy(0.0d0,x_anharm2,
     &       harmfreq2,iWork(iplevel1),
     &       0.0d0,x_anharm2,harmfreq2,iWork(iplevel2),
     &        Work(ipFreqDiffMat+iOrd),l_harm)
      End Do
      Do jOrd = 0,max_nOrd
      dE = Work(ipFreqDiffMat+jOrd)*HarToaJ*1.0d-18
      Do iOrd = 0,max_mOrd
      const2 = const1*exp(-dE/(1.38066d-23*Temperature))
      IntensityMat(iOrd,jOrd) = const2*
     &            Abs(TermMat(iOrd,jOrd))*
     &            (Work(ipFC2+iOrd+(max_mOrd+1)*
     &            (jOrd +(max_nOrd+1)))**2+
     &             Work(ipFC2+iOrd+(max_mOrd+1)*
     &            (jOrd +(max_nOrd+1)*2))**2+
     &             Work(ipFC2+iOrd+(max_mOrd+1)*
     &            (jOrd +(max_nOrd+1)*3))**2)
      If ( .not.OscStr ) Then
      IntensityMat(iOrd,jOrd) = 32.13002d9*
     &             const2*(TermMat(iOrd,jOrd)**2)*
     &             IntensityMat(iOrd,jOrd)
      End If
      End Do
      End Do
      End If
      Call GetMem('level1','Free','Inte',iplevel1,nOsc)
      Call GetMem('level2','Free','Inte',iplevel2,nOsc)
      Call GetMem('FC2','Free','Real',
     &            ipFC2,(max_mOrd+1)*(max_nOrd+1)*4)
      Call GetMem('FreqDiffMat','Free','Real',
     &    ipFreqDiffMat,l_FreqDiffMat)

C!
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(T0)
         Call Unused_integer(max_term)
         Call Unused_integer(nDimTot)
      End If
      End
