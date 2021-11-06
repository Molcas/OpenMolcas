!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1998, Anders Bernhardsson                              *
!               1998, Niclas Forsberg                                  *
!***********************************************************************

!module ForceFieldIntMod

!  Purpose:
!    Calculates transition intensities using the double harmonic
!    approximation.
!
!  Written by:
!    Anders Bernhardsson & Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1998.

!vv private

!contains

subroutine IntForceField(IntensityMat,TermMat,T0,max_term,FC00,C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,r00,m_max,n_max,max_dip, &
                         Trandip,TranDipGrad,harmfreq1,x_anharm1,harmfreq2,x_anharm2,mMat,mInc,mDec,nMat,nInc,nDec,OscStr,nnsiz, &
                         max_mOrd,max_nOrd,nDimTot,nOsc)
!  Purpose:
!    Calculates the intensities of the different transitions between
!    the two surfaces.
!
!  Input:
!    FC           : Real*8 two dimensional array -
!                   Franck-Condon factors.
!    T0           : Real*8 variable - energy difference between
!                   the two states.
!    max_term     : Integer - maximum order of the transition dipole terms.
!    W1,W2        : Real*8 two dimensional arrays - eigenvectors
!                   scaled by the square root of the eigenvalues. Harmonic
!                   approximation.
!    C1,C2        : Real*8 two dimensional arrays - inverses
!                   of W1 and W2.
!    det1,det2    : Real*8 variables - determinants of C1 and C2.
!    r01,r02      : Real*8 arrays - coordinates of the two
!                   oscillators.
!    r00          : Real*8 array - coordinates of intermediate
!                   oscillator.
!    m_max,n_max  : Integer variables - maximum quanta.
!    max_dip      : Integer variable - maximum order of transition dipole.
!    MatEl        : Logical
!    m_plot,
!    n_plot       : Integer array - transitions wanted in output.
!    TermMat      : Real*8 two dimensional array - energies
!                   of transitions.
!
!  Output:
!    IntensityMat : Real*8 two dimensional array - intensities
!                   of the transitions.
!
!  Uses:
!    Constants
!    VibMod
!    TabMod
!    FCMod
!    MatElMod

!use VibMod
!use TabMod
!use FCMod
!use MatElMod
implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
#include "dims.fh"
real*8 IntensityMat(0:max_mOrd,0:max_nOrd)
real*8 TermMat(0:max_mOrd,0:max_nOrd)
real*8 C1(nosc,nosc), C2(nosc,nosc), W1(nosc,nosc), W2(nosc,nosc), C(nosc,nosc), W(nosc,nosc)
real*8 r01(nosc), r02(nosc), r00(nosc)
real*8 TranDip(3)
real*8 TranDipGrad(3,nosc)
real*8 harmfreq1(nosc), harmfreq2(nosc)
real*8 x_anharm1(nosc,nosc), x_anharm2(nosc,nosc)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2)
logical OscStr
integer nvTabDim
#include "WrkSpc.fh"

! Initialize.
call TabDim2_drv(m_max,nosc,nvTabDim)
max_mOrd = nvTabDim-1
call TabDim2_drv(n_max,nosc,nvTabDim)
max_nOrd = nvTabDim-1
call GetMem('FC2','Allo','Real',ipFC2,(max_mOrd+1)*(max_nOrd+1)*4)

call SetUpHarmDip(Work(ipFC2),max_dip,m_max,n_max,mMat,mInc,mDec,nMat,nInc,nDec,C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,r00, &
                  TranDip,TranDipGrad,FC00,nnsiz,max_mOrd,max_nOrd,nOsc)

! Calculate intensities with Boltzmann weighting of hotband intensity.
const1 = (2.0d0/3.0d0)
call GetMem('level1','Allo','Inte',iplevel1,nOsc)
call GetMem('level2','Allo','Inte',iplevel2,nOsc)

if (max_nOrd > max_mOrd) then
  l_FreqDiffMat = max_mOrd+1
  call GetMem('FreqDiffMat','Allo','Real',ipFreqDiffMat,l_FreqDiffMat)

  do iv=1,nOsc
    iWork(iplevel1+iv-1) = mMat(0,iv)
  end do
  do iOrd=0,max_mOrd
    do iv=1,nOsc
      iWork(iplevel2+iv-1) = mMat(iOrd,iv)
    end do
    l_harm = nOsc
    call TransEnergy(0.0d0,x_anharm1,harmfreq1,iWork(iplevel1),0.0d0,x_anharm1,harmfreq1,iWork(iplevel2),Work(ipFreqDiffMat+iOrd), &
                     l_harm)
  end do
  do jOrd=0,max_nOrd
    do iOrd=0,max_mOrd
      dE = Work(ipFreqDiffMat+iOrd)*HarToaJ*1.0d-18
      const2 = const1*exp(-dE/(1.38066d-23*Temperature))
      IntensityMat(iOrd,jOrd) = const2*abs(TermMat(iOrd,jOrd))*(Work(ipFC2+iOrd+(max_mOrd+1)*(jOrd+(max_nOrd+1)))**2+ &
                                Work(ipFC2+iOrd+(max_mOrd+1)*(jOrd+(max_nOrd+1)*2))**2+ &
                                Work(ipFC2+iOrd+(max_mOrd+1)*(jOrd+(max_nOrd+1)*3))**2)
      if (.not. OscStr) then
        IntensityMat(iOrd,jOrd) = 32.13002d9*const2*(TermMat(iOrd,jOrd)**2)*IntensityMat(iOrd,jOrd)
      end if
    end do
  end do
else
  l_FreqDiffMat = max_nOrd+1
  call GetMem('FreqDiffMat','Allo','Real',ipFreqDiffMat,l_FreqDiffMat)
  do iv=1,nOsc
    iWork(iplevel1+iv-1) = nMat(0,iv)
  end do
  do iOrd=0,max_nOrd
    do iv=1,nOsc
      iWork(iplevel2+iv-1) = nMat(iOrd,iv)
    end do
    l_harm = nOsc
    call TransEnergy(0.0d0,x_anharm2,harmfreq2,iWork(iplevel1),0.0d0,x_anharm2,harmfreq2,iWork(iplevel2),Work(ipFreqDiffMat+iOrd), &
                     l_harm)
  end do
  do jOrd=0,max_nOrd
    dE = Work(ipFreqDiffMat+jOrd)*HarToaJ*1.0d-18
    do iOrd=0,max_mOrd
      const2 = const1*exp(-dE/(1.38066d-23*Temperature))
      IntensityMat(iOrd,jOrd) = const2*abs(TermMat(iOrd,jOrd))*(Work(ipFC2+iOrd+(max_mOrd+1)*(jOrd+(max_nOrd+1)))**2+ &
                                Work(ipFC2+iOrd+(max_mOrd+1)*(jOrd+(max_nOrd+1)*2))**2+ &
                                Work(ipFC2+iOrd+(max_mOrd+1)*(jOrd+(max_nOrd+1)*3))**2)
      if (.not. OscStr) then
        IntensityMat(iOrd,jOrd) = 32.13002d9*const2*(TermMat(iOrd,jOrd)**2)*IntensityMat(iOrd,jOrd)
      end if
    end do
  end do
end if
call GetMem('level1','Free','Inte',iplevel1,nOsc)
call GetMem('level2','Free','Inte',iplevel2,nOsc)
call GetMem('FC2','Free','Real',ipFC2,(max_mOrd+1)*(max_nOrd+1)*4)
call GetMem('FreqDiffMat','Free','Real',ipFreqDiffMat,l_FreqDiffMat)

! Avoid unused argument warnings
if (.false.) then
  call Unused_real(T0)
  call Unused_integer(max_term)
  call Unused_integer(nDimTot)
end if

end subroutine IntForceField
