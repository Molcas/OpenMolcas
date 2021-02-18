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
* Copyright (C) 1996, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
c       Module IntensityMod
C!
C!  Contains:
C!    DipMatEl       (Dij,W,L,U,FC00,nMat,D0,D1,D2,D3,D4,max_term)
C!    SetUpDipMat    (DipMat,max_term,ipow,var,dip,trfName,
C!                    use_weight,C1,W1,det1,r01,
C!                    C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,
C!                    max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,
C!                    nInc,mDec,nDec)
C!    Intensity      (IntensityMat,TermMat,T0,harmfreq1,harmfreq2,
C!                    x_anharm1,x_anharm2,max_term,ipow,var,
C!                    Tdip_y,Tdip_z,
C!                    trfName,use_weight,U1,U2,E1,E2,r00,
C!                    C1,W1,det1,r01,C2,W2,det2,r02,m_max,n_max,
C!                    max_dip,MatEl)
C!
C!  Uses:
C!    TabMod
C!    FCMod
C!    MatElMod
C!    OptMod
C!    VibMod
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!
C!-----------------------------------------------------------------------!
C!
cxx       Private
C!
c       Contains




C!-----------------------------------------------------------------------!
C!
      Subroutine Intensity(IntensityMat,TermMat,T0,max_term,ipow,var,
     &    Tdip_x,Tdip_y,Tdip_z,
     &    trfName,U1,U2,E1,E2,C1,W1,det1,r01,
     &    C2,W2,det2,r02,C,W,det0,r00,m_max,n_max,max_dip,
     &    m_plot,n_plot,harmfreq1,harmfreq2,r0,r1,r2,Base,
     &   l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,
     &   nOsc,nDimTot,nPolyTerm,ndata,nvar,MaxNumAt,l_n_plot,l_m_plot)
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
C!    ipow         : Two dimensional integer array - terms of the
C!                   polynomial.
C!    var          : Real*8 two dimensional array - coordinates
C!                   to be used in the fit.
C!    Tdip_y,
C!    Tdip_z       : Real*8 array - transition dipole
C!    trfName      : Character array - transformation associated with each
C!                   internal coordinate.
C!    use_weight   : Logical
C!    U1,U2        : Real*8 two dimensional arrays - eigenvectors
C!                   obtained from matrix element calculations.
C!    E1,E2        : Real*8 arrays - energy values.
C!                   obtained from matrix element calculations.
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
C!    m_plot,
C!    n_plot       : Integer array - transitions wanted in output.
C!
C!  Output:
C!    IntensityMat : Real*8 two dimensional array - intensities
C!                   of the transitions.
C!    TermMat      : Real*8 two dimensional array - energies
C!                   of transitions.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
      Real*8 IntensityMat(0:l_IntensityMat_1,0:l_IntensityMat_2)
      Real*8 TermMat(0:l_TermMat_1,0:l_TermMat_2)
      Integer ipow(nPolyTerm,nvar)
      Real*8 var(ndata,nvar)
      Real*8 Tdip_x(ndata), Tdip_y(ndata),Tdip_z(ndata)
      Character*80 trfName(MaxNumAt)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),
     &  W1(nOsc,nOsc),W2(nOsc,nOsc),C(nOsc,nOsc),W(nOsc,nOsc)
      Real*8 r01(nOsc),r02(nOsc),r00(nOsc)
      Real*8 r1(nOsc),r2(nOsc),r0(nOsc)
      Real*8 U1(nDimTot,nDimTot),U2(nDimTot,nDimTot)
      Real*8 Base (nOsc,nOsc)
      Real*8 E1(nDimTot),E2(nDimTot)
      Real*8 harmfreq1(nOsc),harmfreq2(nOsc)
      Integer m_plot(l_m_plot),n_plot(l_n_plot)

      Integer  nvTabDim
#include "WrkSpc.fh"
C!
C!---- Calculate dimensions given max level of excitation for the
C!     different states.
      Call TabDim2_drv(m_max,nOsc,nvTabDim)
      mTabDim = nvTabDim-1
      Call TabDim2_drv(n_max,nOsc,nvTabDim)
      nTabDim = nvTabDim-1
C! n_max2 = n_max+max_dip
C! max_term = n_max2-n_max
C! nTabDim2 = TabDim2(n_max2,nOsc)-1
      n_max2 = n_max
      Call TabDim2_drv(n_max2,nOsc,nvTabDim)
      nTabDim2 = nvTabDim-1
      max_mOrd = mTabDim
      max_nOrd = nTabDim
      max_nOrd2 = nTabDim
C!
C!---- Set up mMat for L.
      Call GetMem('mMat','Allo','Inte',ipmMat,(mTabDim+1)*nOsc)
      Call GetMem('mInc','Allo','Inte',ipmInc,(mTabDim+1)*nOsc)
      Call GetMem('mDec','Allo','Inte',ipmDec,(mTabDim+1)*nOsc)
c Put dimensions into common block:
      mdim1=mTabDim
      mdim2=nOsc
      Call MakeTab2(
     &   m_max,max_mOrd,max_mInc,mTabDim,iWork(ipmMat),
     &   iWork(ipmInc),iWork(ipmDec),nOsc)
C!
C!---- Set up nMat for U.
      max_nOrd = nTabDim
      Call TabDim2_drv(max(0,n_max-1),nOsc,nvTabDim)
      max_nInc = nvTabDim-1
      Call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim2+1)*nOsc)
      Call GetMem('nInc','Allo','Inte',ipnInc,(nTabDim2+1)*nOsc)
      Call GetMem('nDec','Allo','Inte',ipnDec,(nTabDim2+1)*nOsc)
c Put dimensions into common block:
      ndim1=nTabDim2
      ndim2=nOsc
c Use nnsiz for the time being, to transfer dim to called routines.
      nnsiz=ntabdim2
      Call MakeTab2(
     &     n_max2,max_nOrd2,max_nInc2,nTabDim2,iWork(ipnMat),
     &     iWork(ipnInc),iWork(ipnDec),nOsc)
C!
C!---- Either use the eigenvectors obtained from the variational
C!     method or use the simpler harmonic approximation.
      nDimTot = 2*max_mOrd+2
      n=ndimtot-1
      Call GetMem('FC2','Allo','Real',ipFC2,nDimTot*nDimTot*3)
c       FC2=0.0d0
      call dcopy_(nDimtot*nDimTot*3,[0.0d0],0,Work(ipFC2),1)
      Call GetMem('DipMat','Allo','Real',ipDipMat,nDimTot*nDimTot*3)
c       DipMat = 0.0d0
      call dcopy_(nDimtot*nDimTot*3,[0.0d0],0,Work(ipDipMat),1)
      Call SetUpDipMat(Work(ipDipMat),max_dip,ipow,var,
     &     Tdip_x,trfName,C1,W1,det1,r01,
     &     C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,
     &     max_mInc,max_nInc,max_nInc,iWork(ipmMat),
     &     iWork(ipnMat),iWork(ipmInc),
     &     iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &     det0,r0,r1,r2,base,nnsiz,nOsc,
     &   nDimTot,nPolyTerm,ndata,nvar,MaxNumAt)
      Call SetUpDipMat(Work(ipDipMat+nDimTot*nDimTot),
     &  max_dip,ipow,var,
     &    Tdip_y,trfName,C1,W1,det1,r01,
     &    C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,
     &    max_mInc,max_nInc,max_nInc,iWork(ipmMat),
     &    iWork(ipnMat),iWork(ipmInc),
     &    iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &    det0,r0,r1,r2,base,nnsiz,nOsc,
     &   nDimTot,nPolyTerm,ndata,nvar,MaxNumAt)
      m_plot_max = l_m_plot
      n_plot_max = l_n_plot
      max_mQuanta = m_plot(1)
      If ( m_plot_max.gt.1 ) Then
      Do i = 2,m_plot_max
      If ( m_plot(i).gt.max_mQuanta ) Then
      max_mQuanta = m_plot(i)
      End If
      End Do
      End If
      max_nQuanta = n_plot(1)
      If ( n_plot_max.gt.1 ) Then
      Do i = 2,n_plot_max
      If ( n_plot(i).gt.max_nQuanta ) Then
      max_nQuanta = n_plot(i)
      End If
      End Do
      End If
C!  max_mOrd = nDimTot-1 ; max_nOrd = nDimTot-1
C!  If ( max_mQuanta.eq.0 ) max_mOrd = 0
C!  If ( max_nQuanta.eq.0 ) max_nOrd = 0
      Do k = 1,2
      Do m = 1,ndimtot
      Do n = 1,ndimtot
      sum = 0.0d0
      Do j = 1,nDimTot
      Do i = 1,nDimTot
      sum = sum+
     &          Work(ipDipMat+i-1+nDimTot*((j-1)+nDimTot*(k-1)))*
     &          U1(i,m)*U2(j,n)
      End Do
      End Do
      Work(ipFC2+m-1+nDimTot*((n-1)+nDimTot*(k-1))) = sum
      End Do
      End Do
      End Do
      Call GetMem('DipMat','Free','Real',ipDipMat,nDimTot*nDimTot*3)

C!
C!---- Calculate frequency differences to be used in the boltzmann
C!     weighting of the different transitions.
C! level1 = mMat(0,:)
C! Do iOrd = 0,max_mOrd
C!    level2 = mMat(iOrd,:)
C!    FreqDiffMat(iOrd) = TransEnergy(0.0d0,x_anharm1,harmfreq1,level1, &
C!                                    0.0d0,x_anharm1,harmfreq1,level2)
C! End Do
C!
C!---- Calculate intensities.
      const1 = (2.0d0/3.0d0)*32.13002d9
c       Real*8 IntensityMat(0:l_IntensityMat_1,0:l_IntensityMat_2)
      call dcopy_((l_IntensityMat_1+1)*(l_IntensityMat_2+1),
     &  [0.0d0],0,IntensityMat,1)
c       IntensityMat = 0.0d0
      Do jOrd = 0,ndimtot-1
      Do iOrd = 0,ndimtot-1
C!      dE = FreqDiffMat(iOrd)*HarToaJ*1.0d-18
C!      const2 = const1*exp(-dE/(1.38066d-23*Temperature))
C!      IntensityMat(iOrd,jOrd) = const2*(TermMat(iOrd,jOrd)**3)* &
      IntensityMat(iOrd,jOrd) = const1*
     &       (TermMat(iOrd,jOrd)**3)*
     &                  (Work(ipFC2+iOrd+nDimTot*(jOrd))**2+
     &      Work(ipFC2+iOrd+nDimTot*(jOrd+nDimTot))**2+
     &     Work(ipFC2+iOrd+nDimTot*(jOrd+nDimTot*2))**2)
      End Do
      End Do

C!
      Call GetMem('FC2','Free','Real',ipFC2,nDimTot*nDimTot*3)

      Call GetMem('mMat','Free','Inte',ipmMat,(mTabDim+1)*nOsc)
      Call GetMem('mInc','Free','Inte',ipmInc,(mTabDim+1)*nOsc)
      Call GetMem('mDec','Free','Inte',ipmDec,(mTabDim+1)*nOsc)
      Call GetMem('nMat','Free','Inte',ipnMat,(nTabDim2+1)*nOsc)
      Call GetMem('nInc','Free','Inte',ipnInc,(nTabDim2+1)*nOsc)
      Call GetMem('nDec','Free','Inte',ipnDec,(nTabDim2+1)*nOsc)
C!
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(T0)
         Call Unused_integer(max_term)
         Call Unused_real_array(Tdip_z)
         Call Unused_real_array(E1)
         Call Unused_real_array(E2)
         Call Unused_real_array(C)
         Call Unused_real_array(W)
         Call Unused_real_array(r00)
         Call Unused_real_array(harmfreq1)
         Call Unused_real_array(harmfreq2)
      End If
      End
C!

C!-----------------------------------------------------------------------!


C!-----------------------------------------------------------------------!
C!
      Subroutine Intensity2(IntensityMat,TermMat,T0,max_term,
     &  U1,U2,E1,E2,C1,W1,det1,r01,
     &  C2,W2,det2,r02,C,W,det0,r00,m_max,n_max,max_dip,
     &  m_plot,n_plot,TranDip,TranDipGrad,
     &  harmfreq1,x_anharm1,harmfreq2,x_anharm2,
     &  r0,r1,r2,Base,
     &  l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,
     &  nOsc,nDimTot,l_n_plot,l_m_plot)
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
C!    ipow         : Two dimensional integer array - terms of the
C!                   polynomial.
C!    var          : Real*8 two dimensional array - coordinates
C!                   to be used in the fit.
C!    Tdip_y,
C!    Tdip_z       : Real*8 array - transition dipole
C!    trfName      : Character array - transformation associated with each
C!                   internal coordinate.
C!    use_weight   : Logical
C!    U1,U2        : Real*8 two dimensional arrays - eigenvectors
C!                   obtained from matrix element calculations.
C!    E1,E2        : Real*8 arrays - energy values.
C!                   obtained from matrix element calculations.
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
C!    m_plot,
C!    n_plot       : Integer array - transitions wanted in output.
C!
C!  Output:
C!    IntensityMat : Real*8 two dimensional array - intensities
C!                   of the transitions.
C!    TermMat      : Real*8 two dimensional array - energies
C!                   of transitions.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
      Real*8 IntensityMat(0:l_IntensityMat_1,0:l_IntensityMat_2)
      Real*8 TermMat(0:l_TermMat_1,0:l_TermMat_2)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),
     &  W1(nOsc,nOsc),W2(nOsc,nOsc),C(nOsc,nOsc),W(nOsc,nOsc)
      Real*8 r01(nOsc),r02(nOsc),r00(nOsc)
      Real*8 r1(nOsc),r2(nOsc),r0(nOsc)
      Real*8 U1(nDimTot,nDimTot),U2(nDimTot,nDimTot)
      Real*8 TranDip(3)
      Real*8 TranDipGrad(3,nOsc)
      Real*8 Base (nOsc,nOsc)
      Real*8 E1(nDimTot),E2(nDimTot)
      Real*8 harmfreq1(nOsc),harmfreq2(nOsc)
      Real*8 x_anharm1(nOsc,nOsc),x_anharm2(nOsc,nOsc)
      Integer m_plot(l_m_plot),n_plot(l_n_plot)

      Integer  nvTabDim
#include "WrkSpc.fh"
C!
C!---- Calculate dimensions given max level of excitation for the
C!     different states.
      Call TabDim2_drv(m_max,nOsc,nvTabDim)
      mTabDim = nvTabDim-1
      Call TabDim2_drv(n_max,nOsc,nvTabDim)
      nTabDim = nvTabDim-1
      n_max2 = n_max
      Call TabDim2_drv(n_max2,nOsc,nvTabDim)
      nTabDim2 = nvTabDim-1
      max_mOrd = mTabDim
      max_nOrd = nTabDim
      max_nOrd2 = nTabDim
C!
C!---- Set up mMat for L.
      Call GetMem('mMat','Allo','Inte',ipmMat,(mTabDim+1)*nOsc)
      Call GetMem('mInc','Allo','Inte',ipmInc,(mTabDim+1)*nOsc)
      Call GetMem('mDec','Allo','Inte',ipmDec,(mTabDim+1)*nOsc)
c Put dimensions into common block:
      mdim1=mTabDim
      mdim2=nOsc
      Call MakeTab2(
     &   m_max,max_mOrd,max_mInc,mTabDim,iWork(ipmMat),
     &   iWork(ipmInc),iWork(ipmDec),nOsc)
C!
C!---- Set up nMat for U.
      max_nOrd = nTabDim
      Call TabDim2_drv(max(0,n_max-1),nOsc,nvTabDim)
      max_nInc = nvTabDim-1
      Call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim2+1)*nOsc)
      Call GetMem('nInc','Allo','Inte',ipnInc,(nTabDim2+1)*nOsc)
      Call GetMem('nDec','Allo','Inte',ipnDec,(nTabDim2+1)*nOsc)
c Put dimensions into common block:
      ndim1=nTabDim2
      ndim2=nOsc
      nnsiz=ntabdim2
      Call MakeTab2(
     &  n_max2,max_nOrd2,max_nInc2,nTabDim2,iWork(ipnMat),
     &  iWork(ipnInc),iWork(ipnDec),nOsc)
C!
C!---- Either use the eigenvectors obtained from the variational
C!     method or use the simpler harmonic approximation.
      nDimTot = max_mOrd+1
      n = nDimTot-1
      Call GetMem('FC2','Allo','Real',ipFC2,nDimTot*nDimTot*3)
c       FC2 = 0.0d0
      call dcopy_(nDimtot*nDimTot*3,[0.0d0],0,Work(ipFC2),1)
      Call GetMem('DipMat','Allo','Real',ipDipMat,nDimTot*nDimTot*3)
c       DipMat = 0.0d0
      call dcopy_(nDimtot*nDimTot*3,[0.0d0],0,Work(ipDipMat),1)
      Call SetUpDipMat2(Work(ipDipMat),max_dip,
     &       C1,W1,det1,r01,C2,W2,det2,r02,
     &       max_mOrd,max_nOrd,max_nOrd,
     &       max_mInc,max_nInc,max_nInc,
     &       iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),
     &       iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &       det0,r0,r1,r2,Base,TranDip(1),TranDipGrad(1,1),
     &       nnsiz,nOsc,nDimTot)
      Call SetUpDipMat2(Work(ipDipMat+nDimTot*nDimTot),max_dip,
     &       C1,W1,det1,r01,C2,W2,det2,r02,
     &       max_mOrd,max_nOrd,max_nOrd,
     &       max_mInc,max_nInc,max_nInc,
     &       iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),
     &       iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &       det0,r0,r1,r2,Base,TranDip(2),TranDipGrad(2,1),
     &       nnsiz,nOsc,nDimTot)
      Call SetUpDipMat2(Work(ipDipMat+0+nDimTot*nDimTot*2),max_dip,
     &       C1,W1,det1,r01,C2,W2,det2,r02,
     &       max_mOrd,max_nOrd,max_nOrd,
     &       max_mInc,max_nInc,max_nInc,
     &       iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),
     &       iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &       det0,r0,r1,r2,Base,TranDip(3),TranDipGrad(3,1),
     &       nnsiz,nOsc,nDimTot)
      m_plot_max = l_m_plot
      n_plot_max = l_n_plot
      max_mQuanta = m_plot(1)
      If ( m_plot_max.gt.1 ) Then
      Do i = 2,m_plot_max
      If ( m_plot(i).gt.max_mQuanta ) Then
      max_mQuanta = m_plot(i)
      End If
      End Do
      End If
      max_nQuanta = n_plot(1)
      If ( n_plot_max.gt.1 ) Then
      Do i = 2,n_plot_max
      If ( n_plot(i).gt.max_nQuanta ) Then
      max_nQuanta = n_plot(i)
      End If
      End Do
      End If
C! max_mOrd = nDimTot-1 ; max_nOrd = nDimTot-1
C! If ( max_mQuanta.eq.0 ) max_mOrd = 0
C! If ( max_nQuanta.eq.0 ) max_nOrd = 0
      Do k = 1,3
      Do m = 1,nDimTot
      Do n = 1,nDimTot
      sum = 0.0d0
      Do j = 1,nDimTot
      Do i = 1,nDimTot
      sum = sum+
     &                Work(ipDipMat+(i-1)+nDimTot*((j-1)+
     &                nDimTot*(k-1)))*U1(i,m)*U2(j,n)
      End Do
      End Do
      Work(ipFC2+m-1+nDimTot*((n-1)+nDimTot*(k-1))) = sum
      End Do
      End Do
      End Do
      Call GetMem('DipMat','Free','Real',ipDipMat,nDimTot*nDimTot*3)
C!
C!---- Calculate frequency differences to be used in the boltzmann
C!     weighting of the different transitions.
C! level1 = mMat(0,:)
C! Do iOrd = 0,max_mOrd
C!    level2 = mMat(iOrd,:)
C!    FreqDiffMat(iOrd) = TransEnergy(0.0d0,x_anharm1,harmfreq1,level1, &
C!                                    0.0d0,x_anharm1,harmfreq1,level2)
C! End Do
C!
C!---- Calculate intensities.
      const1 = (2.0d0/3.0d0)*32.13002d9
      call dcopy_((l_IntensityMat_1+1)*(l_IntensityMat_2+1),
     &  [0.0d0],0,IntensityMat,1)
c       IntensityMat = 0.0d0
      Do jOrd = 0,nDimTot-1
      Do iOrd = 0,nDimTot-1
C!      dE = FreqDiffMat(iOrd)*HarToaJ*1.0d-18
C!      const2 = const1*exp(-dE/(1.38066d-23*Temperature))
C!      IntensityMat(iOrd,jOrd) = const2*(TermMat(iOrd,jOrd)**3)* &
      IntensityMat(iOrd,jOrd) = const1*
     &       (TermMat(iOrd,jOrd)**3)*
     &                  (Work(ipFC2+iOrd+nDimTot*(jOrd))**2+
     &     Work(ipFC2+iOrd+nDimTot*(jOrd+nDimTot))**2+
     &     Work(ipFC2+iOrd+nDimTot*(jOrd+nDimTot*2))**2)
      End Do
      End Do

C!
      Call GetMem('FC2','Free','Real',ipFC2,nDimTot*nDimTot*3)
      Call GetMem('mMat','Free','Inte',ipmMat,(mTabDim+1)*nOsc)
      Call GetMem('mInc','Free','Inte',ipmInc,(mTabDim+1)*nOsc)
      Call GetMem('mDec','Free','Inte',ipmDec,(mTabDim+1)*nOsc)
      Call GetMem('nMat','Free','Inte',ipnMat,(nTabDim2+1)*nOsc)
      Call GetMem('nInc','Free','Inte',ipnInc,(nTabDim2+1)*nOsc)
      Call GetMem('nDec','Free','Inte',ipnDec,(nTabDim2+1)*nOsc)
C!
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(T0)
         Call Unused_integer(max_term)
         Call Unused_real_array(E1)
         Call Unused_real_array(E2)
         Call Unused_real_array(C)
         Call Unused_real_array(W)
         Call Unused_real_array(r00)
         Call Unused_real_array(harmfreq1)
         Call Unused_real_array(x_anharm1)
         Call Unused_real_array(harmfreq2)
         Call Unused_real_array(x_anharm2)
      End If
      End
