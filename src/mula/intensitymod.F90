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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

!module IntensityMod

!  Contains:
!    DipMatEl       (Dij,W,L,U,FC00,nMat,D0,D1,D2,D3,D4,max_term)
!    SetUpDipMat    (DipMat,max_term,ipow,var,dip,trfName,
!                    use_weight,C1,W1,det1,r01,
!                    C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,
!                    max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,
!                    nInc,mDec,nDec)
!    Intensity      (IntensityMat,TermMat,T0,harmfreq1,harmfreq2,
!                    x_anharm1,x_anharm2,max_term,ipow,var,
!                    Tdip_y,Tdip_z,
!                    trfName,use_weight,U1,U2,E1,E2,r00,
!                    C1,W1,det1,r01,C2,W2,det2,r02,m_max,n_max,
!                    max_dip,MatEl)
!
!  Uses:
!    TabMod
!    FCMod
!    MatElMod
!    OptMod
!    VibMod
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

!xx private

!contains

subroutine Intensity(IntensityMat,TermMat,ipow,var,Tdip_x,Tdip_y,trfName,U1,U2,C1,W1,det1,r01,C2,W2,det2,r02,det0,m_max,n_max, &
                     max_dip,m_plot,n_plot,r0,r1,r2,Base,l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc,nDimTot, &
                     nPolyTerm,ndata,nvar,MaxNumAt,l_n_plot,l_m_plot)
!  Purpose:
!    Calculates the intensities of the different transitions between
!    the two surfaces.
!
!  Input:
!    FC           : Real*8 two dimensional array -
!                   Franck-Condon factors.
!    ipow         : Two dimensional integer array - terms of the
!                   polynomial.
!    var          : Real*8 two dimensional array - coordinates
!                   to be used in the fit.
!    Tdip_x,
!    Tdip_y       : Real*8 array - transition dipole
!    trfName      : Character array - transformation associated with each
!                   internal coordinate.
!    use_weight   : Logical
!    U1,U2        : Real*8 two dimensional arrays - eigenvectors
!                   obtained from matrix element calculations.
!    W1,W2        : Real*8 two dimensional arrays - eigenvectors
!                   scaled by the square root of the eigenvalues. Harmonic
!                   approximation.
!    C1,C2        : Real*8 two dimensional arrays - inverses
!                   of W1 and W2.
!    det1,det2    : Real*8 variables - determinants of C1 and C2.
!    r01,r02      : Real*8 arrays - coordinates of the two
!                   oscillators.
!    m_max,n_max  : Integer variables - maximum quanta.
!    max_dip      : Integer variable - maximum order of transition dipole.
!    m_plot,
!    n_plot       : Integer array - transitions wanted in output.
!
!  Output:
!    IntensityMat : Real*8 two dimensional array - intensities
!                   of the transitions.
!    TermMat      : Real*8 two dimensional array - energies
!                   of transitions.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three
use Definitions, only: wp

implicit real*8(a-h,o-z)
real*8 IntensityMat(0:l_IntensityMat_1,0:l_IntensityMat_2)
real*8 TermMat(0:l_TermMat_1,0:l_TermMat_2)
integer ipow(nPolyTerm,nvar)
real*8 var(ndata,nvar)
real*8 Tdip_x(ndata), Tdip_y(ndata)
character*80 trfName(MaxNumAt)
real*8 C1(nOsc,nOsc), C2(nOsc,nOsc), W1(nOsc,nOsc), W2(nOsc,nOsc)
real*8 r01(nOsc), r02(nOsc)
real*8 r1(nOsc), r2(nOsc), r0(nOsc)
real*8 U1(nDimTot,nDimTot), U2(nDimTot,nDimTot)
real*8 Base(nOsc,nOsc)
integer m_plot(l_m_plot), n_plot(l_n_plot)
integer nvTabDim
integer, allocatable :: mDec(:,:), mInc(:,:), mMat(:,:), nDec(:,:), nInc(:,:), nMat(:,:)
real*8, allocatable :: DipMat(:,:,:), FC2(:,:,:)

! Calculate dimensions given max level of excitation for the different states.
call TabDim(m_max,nOsc,nvTabDim)
mTabDim = nvTabDim-1
call TabDim(n_max,nOsc,nvTabDim)
nTabDim = nvTabDim-1
!n_max2 = n_max+max_dip
!max_term = n_max2-n_max
!nTabDim2 = TabDim2(n_max2,nOsc)-1
n_max2 = n_max
call TabDim(n_max2,nOsc,nvTabDim)
nTabDim2 = nvTabDim-1
max_mOrd = mTabDim
max_nOrd = nTabDim
max_nOrd2 = nTabDim

! Set up mMat for L.
call mma_allocate(mMat,[0,mTabDim],[1,nOsc],label='mMat')
call mma_allocate(mInc,[0,mTabDim],[1,nOsc],label='mInc')
call mma_allocate(mDec,[0,mTabDim],[1,nOsc],label='mDec')
! Put dimensions into module
mdim1 = mTabDim
mdim2 = nOsc
call MakeTab2(m_max,max_mOrd,max_mInc,mTabDim,mMat,mInc,mDec,nOsc)

! Set up nMat for U.
max_nOrd = nTabDim
call TabDim(max(0,n_max-1),nOsc,nvTabDim)
max_nInc = nvTabDim-1
call mma_allocate(nMat,[0,nTabDim],[1,nOsc],label='nMat')
call mma_allocate(nInc,[0,nTabDim],[1,nOsc],label='nInc')
call mma_allocate(nDec,[0,nTabDim],[1,nOsc],label='nDec')
! Put dimensions into module
ndim1 = nTabDim2
ndim2 = nOsc
call MakeTab2(n_max2,max_nOrd2,max_nInc2,nTabDim2,nMat,nInc,nDec,nOsc)

! Either use the eigenvectors obtained from the variational
! method or use the simpler harmonic approximation.
nDimTot = 2*max_mOrd+2
call mma_allocate(FC2,nDimTot,nDimTot,3)
FC2(:,:,:) = Zero
call mma_allocate(DipMat,nDimTot,nDimTot,2,label='DipMat')
DipMat(:,:,:) = Zero
call SetUpDipMat(DipMat(:,:,1),max_dip,ipow,var,Tdip_x,trfName,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,max_mInc, &
                 max_nInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,det0,r0,r1,r2,base,nOsc,nDimTot-1,nPolyTerm,ndata,nvar,MaxNumAt)
call SetUpDipMat(DipMat(:,:,2),max_dip,ipow,var,Tdip_y,trfName,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,max_mInc, &
                 max_nInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,det0,r0,r1,r2,base,nOsc,nDimTot-1,nPolyTerm,ndata,nvar,MaxNumAt)
m_plot_max = l_m_plot
n_plot_max = l_n_plot
max_mQuanta = m_plot(1)
if (m_plot_max > 1) then
  do i=2,m_plot_max
    if (m_plot(i) > max_mQuanta) then
      max_mQuanta = m_plot(i)
    end if
  end do
end if
max_nQuanta = n_plot(1)
if (n_plot_max > 1) then
  do i=2,n_plot_max
    if (n_plot(i) > max_nQuanta) then
      max_nQuanta = n_plot(i)
    end if
  end do
end if
!max_mOrd = nDimTot-1
!max_nOrd = nDimTot-1
!if (max_mQuanta == 0) max_mOrd = 0
!if (max_nQuanta == 0) max_nOrd = 0
do k=1,2
  do m=1,nDimTot
    do n=1,nDimTot
      sum = Zero
      do j=1,nDimTot
        do i=1,nDimTot
          sum = sum+DipMat(i,j,k)*U1(i,m)*U2(j,n)
        end do
      end do
      FC2(m,n,k) = sum
    end do
  end do
end do
call mma_deallocate(DipMat)

! Calculate frequency differences to be used in the boltzmann
! weighting of the different transitions.
!level1 = mMat(0,:)
!do iOrd=0,max_mOrd
!  level2 = mMat(iOrd,:)
!  FreqDiffMat(iOrd) = TransEnergy(Zero,x_anharm1,harmfreq1,level1,Zero,x_anharm1,harmfreq1,level2)
!end do

! Calculate intensities.
! where does this number come from?
const1 = (Two/Three)*32.13002e9_wp
IntensityMat(:,:) = Zero
do jOrd=0,nDimTot-1
  do iOrd=0,nDimTot-1
    !dE = FreqDiffMat(iOrd)*auTokJ*1.0e-39_wp
    !const2 = const1*exp(-dE/(kBoltzmann*Temperature))
    !IntensityMat(iOrd,jOrd) = const2*(TermMat(iOrd,jOrd)**3)* &
    IntensityMat(iOrd,jOrd) = const1*(TermMat(iOrd,jOrd)**3)* &
                              (FC2(iOrd+1,jOrd+1,1)**2+FC2(iOrd+1,jOrd+1,2)**2+FC2(iOrd+1,jOrd+1,3)**2)
  end do
end do

call mma_deallocate(FC2)
call mma_deallocate(mMat)
call mma_deallocate(mInc)
call mma_deallocate(mDec)
call mma_deallocate(nMat)
call mma_deallocate(nInc)
call mma_deallocate(nDec)

end subroutine Intensity
!####
subroutine Intensity2(IntensityMat,TermMat,U1,U2,C1,W1,det1,r01,C2,W2,det2,r02,det0,m_max,n_max,max_dip,m_plot,n_plot,TranDip, &
                      TranDipGrad,Base,l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc,nDimTot,l_n_plot,l_m_plot)

!  Purpose:
!    Calculates the intensities of the different transitions between
!    the two surfaces.
!
!  Input:
!    FC           : Real*8 two dimensional array -
!                   Franck-Condon factors.
!    ipow         : Two dimensional integer array - terms of the
!                   polynomial.
!    var          : Real*8 two dimensional array - coordinates
!                   to be used in the fit.
!    trfName      : Character array - transformation associated with each
!                   internal coordinate.
!    use_weight   : Logical
!    U1,U2        : Real*8 two dimensional arrays - eigenvectors
!                   obtained from matrix element calculations.
!    W1,W2        : Real*8 two dimensional arrays - eigenvectors
!                   scaled by the square root of the eigenvalues. Harmonic
!                   approximation.
!    C1,C2        : Real*8 two dimensional arrays - inverses
!                   of W1 and W2.
!    det1,det2    : Real*8 variables - determinants of C1 and C2.
!    r01,r02      : Real*8 arrays - coordinates of the two
!                   oscillators.
!    m_max,n_max  : Integer variables - maximum quanta.
!    max_dip      : Integer variable - maximum order of transition dipole.
!    m_plot,
!    n_plot       : Integer array - transitions wanted in output.
!
!  Output:
!    IntensityMat : Real*8 two dimensional array - intensities
!                   of the transitions.
!    TermMat      : Real*8 two dimensional array - energies
!                   of transitions.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three
use Definitions, only: wp

implicit real*8(a-h,o-z)
real*8 IntensityMat(0:l_IntensityMat_1,0:l_IntensityMat_2)
real*8 TermMat(0:l_TermMat_1,0:l_TermMat_2)
real*8 C1(nOsc,nOsc), C2(nOsc,nOsc), W1(nOsc,nOsc), W2(nOsc,nOsc)
real*8 r01(nOsc), r02(nOsc)
real*8 U1(nDimTot,nDimTot), U2(nDimTot,nDimTot)
real*8 TranDip(3)
real*8 TranDipGrad(3,nOsc)
real*8 Base(nOsc,nOsc)
integer m_plot(l_m_plot), n_plot(l_n_plot)
integer nvTabDim
integer, allocatable :: mDec(:,:), mInc(:,:), mMat(:,:), nDec(:,:), nInc(:,:), nMat(:,:)
real*8, allocatable :: DipMat(:,:,:), FC2(:,:,:)

! Calculate dimensions given max level of excitation for the different states.
call TabDim(m_max,nOsc,nvTabDim)
mTabDim = nvTabDim-1
call TabDim(n_max,nOsc,nvTabDim)
nTabDim = nvTabDim-1
n_max2 = n_max
call TabDim(n_max2,nOsc,nvTabDim)
nTabDim2 = nvTabDim-1
max_mOrd = mTabDim
max_nOrd = nTabDim
max_nOrd2 = nTabDim

! Set up mMat for L.
call mma_allocate(mMat,[0,mTabDim],[1,nOsc],label='mMat')
call mma_allocate(mInc,[0,mTabDim],[1,nOsc],label='mInc')
call mma_allocate(mDec,[0,mTabDim],[1,nOsc],label='mDec')
! Put dimensions into module
mdim1 = mTabDim
mdim2 = nOsc
call MakeTab2(m_max,max_mOrd,max_mInc,mTabDim,mMat,mInc,mDec,nOsc)

! Set up nMat for U.
max_nOrd = nTabDim
call TabDim(max(0,n_max-1),nOsc,nvTabDim)
max_nInc = nvTabDim-1
call mma_allocate(nMat,[0,nTabDim2],[1,nOsc],label='nMat')
call mma_allocate(nInc,[0,nTabDim2],[1,nOsc],label='nInc')
call mma_allocate(nDec,[0,nTabDim2],[1,nOsc],label='nDec')
! Put dimensions into module
ndim1 = nTabDim2
ndim2 = nOsc
call MakeTab2(n_max2,max_nOrd2,max_nInc2,nTabDim2,nMat,nInc,nDec,nOsc)

! Either use the eigenvectors obtained from the variational
! method or use the simpler harmonic approximation.
nDimTot = max_mOrd+1
call mma_allocate(FC2,nDimTot,nDimTot,3)
FC2(:,:,:) = Zero
call mma_allocate(DipMat,nDimTot,nDimTot,3,label='DipMat')
DipMat(:,:,:) = Zero
call SetUpDipMat2(DipMat(:,:,1),max_dip,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,max_nInc,mMat, &
                  nMat,mInc,nInc,mDec,nDec,det0,Base,TranDip(1),TranDipGrad(1,1),nOsc,nDimTot-1)
call SetUpDipMat2(DipMat(:,:,2),max_dip,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,max_nInc,mMat, &
                  nMat,mInc,nInc,mDec,nDec,det0,Base,TranDip(2),TranDipGrad(2,1),nOsc,nDimTot-1)
call SetUpDipMat2(DipMat(:,:,3),max_dip,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,max_nInc,mMat, &
                  nMat,mInc,nInc,mDec,nDec,det0,Base,TranDip(3),TranDipGrad(3,1),nOsc,nDimTot-1)
m_plot_max = l_m_plot
n_plot_max = l_n_plot
max_mQuanta = m_plot(1)
if (m_plot_max > 1) then
  do i=2,m_plot_max
    if (m_plot(i) > max_mQuanta) then
      max_mQuanta = m_plot(i)
    end if
  end do
end if
max_nQuanta = n_plot(1)
if (n_plot_max > 1) then
  do i=2,n_plot_max
    if (n_plot(i) > max_nQuanta) then
      max_nQuanta = n_plot(i)
    end if
  end do
end if
!max_mOrd = nDimTot-1
!max_nOrd = nDimTot-1
!if (max_mQuanta == 0) max_mOrd = 0
!if (max_nQuanta == 0) max_nOrd = 0
do k=1,3
  do m=1,nDimTot
    do n=1,nDimTot
      sum = Zero
      do j=1,nDimTot
        do i=1,nDimTot
          sum = sum+DipMat(i,j,k)*U1(i,m)*U2(j,n)
        end do
      end do
      FC2(m,n,k) = sum
    end do
  end do
end do
call mma_deallocate(DipMat)

! Calculate frequency differences to be used in the boltzmann
! weighting of the different transitions.
!level1 = mMat(0,:)
!do iOrd=0,max_mOrd
!  level2 = mMat(iOrd,:)
!  FreqDiffMat(iOrd) = TransEnergy(Zero,x_anharm1,harmfreq1,level1,Zero,x_anharm1,harmfreq1,level2)
!end do

! Calculate intensities.
! where does this number come from?
const1 = (Two/Three)*32.13002e9_wp
IntensityMat(:,:) = Zero
do jOrd=0,nDimTot-1
  do iOrd=0,nDimTot-1
    !dE = FreqDiffMat(iOrd)*auTokJ*1.0e-39_wp
    !const2 = const1*exp(-dE/(kBoltzmann*Temperature))
    !IntensityMat(iOrd,jOrd) = const2*(TermMat(iOrd,jOrd)**3)* &
    IntensityMat(iOrd,jOrd) = const1*(TermMat(iOrd,jOrd)**3)* &
                              (FC2(iOrd+1,jOrd+1,1)**2+FC2(iOrd+1,jOrd+1,2)**2+FC2(iOrd+1,jOrd+1,3)**2)
  end do
end do

call mma_deallocate(FC2)
call mma_deallocate(mMat)
call mma_deallocate(mInc)
call mma_deallocate(mDec)
call mma_deallocate(nMat)
call mma_deallocate(nInc)
call mma_deallocate(nDec)

end subroutine Intensity2
