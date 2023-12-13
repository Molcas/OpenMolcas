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
!    SetUpDipMat    (DipMat,max_term,ipow,var,dip,trfName,use_weight,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,
!                    max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec)
!    Intensity      (IntensityMat,TermMat,T0,harmfreq1,harmfreq2,x_anharm1,x_anharm2,max_term,ipow,var,Tdip_y,Tdip_z,trfName,
!                    use_weight,U1,U2,E1,E2,r00,C1,W1,det1,r01,C2,W2,det2,r02,m_max,n_max,max_dip,MatEl)
!    Intensity2     (IntensityMat,TermMat,U1,U2,C1,W1,det1,r01,C2,W2,det2,r02,det0,m_max,n_max,max_dip,m_plot,n_plot,TranDip,
!                    TranDipGrad,Base)
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

!contains

subroutine Intensity(IntensityMat,TermMat,ipow,var,Tdip_x,Tdip_y,trfName,U1,U2,C1,W1,det1,r01,C2,W2,det2,r02,det0,m_max,n_max, &
                     max_dip,m_plot,n_plot,r0,r1,r2,Base,l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc,nDimTot, &
                     nPolyTerm,ndata,nvar,MaxNumAt,l_n_plot,l_m_plot)
!  Purpose:
!    Calculates the intensities of the different transitions between the two surfaces.
!
!  Input:
!    FC           : Real two dimensional array - Franck-Condon factors.
!    ipow         : Two dimensional integer array - terms of the polynomial.
!    var          : Real two dimensional array - coordinates to be used in the fit.
!    Tdip_x,
!    Tdip_y       : Real array - transition dipole
!    trfName      : Character array - transformation associated with each internal coordinate.
!    use_weight   : Logical
!    U1,U2        : Real two dimensional arrays - eigenvectors obtained from matrix element calculations.
!    W1,W2        : Real two dimensional arrays - eigenvectors scaled by the square root of the eigenvalues. Harmonic approximation.
!    C1,C2        : Real two dimensional arrays - inverses of W1 and W2.
!    det1,det2    : Real variables - determinants of C1 and C2.
!    r01,r02      : Real arrays - coordinates of the two oscillators.
!    m_max,n_max  : Integer variables - maximum quanta.
!    max_dip      : Integer variable - maximum order of transition dipole.
!    m_plot,
!    n_plot       : Integer array - transitions wanted in output.
!
!  Output:
!    IntensityMat : Real two dimensional array - intensities of the transitions.
!    TermMat      : Real two dimensional array - energies of transitions.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPolyTerm, nvar, ipow(nPolyTerm,nvar), m_max, n_max, max_dip, l_n_plot, l_m_plot, &
                                 m_plot(l_m_plot), n_plot(l_n_plot), l_IntensityMat_1, l_IntensityMat_2, l_TermMat_1, l_TermMat_2, &
                                 nOsc, ndata, MaxNumAt
real(kind=wp), intent(out) :: IntensityMat(0:l_IntensityMat_1,0:l_IntensityMat_2), det0, r0(nOsc), r1(nOsc), r2(nOsc)
integer(kind=iwp), intent(inout) :: nDimTot
real(kind=wp), intent(in) :: TermMat(0:l_TermMat_1,0:l_TermMat_2), var(ndata,nvar), Tdip_x(ndata), Tdip_y(ndata), &
                             U1(nDimTot,nDimTot), U2(nDimTot,nDimTot), C1(nOsc,nOsc), W1(nOsc,nOsc), det1, r01(nOsc), &
                             C2(nOsc,nOsc), W2(nOsc,nOsc), det2, r02(nOsc), Base(nOsc,nOsc)
character(len=80), intent(in) :: trfName(MaxNumAt)
integer(kind=iwp) :: i, iOrd, j, jOrd, k, m, m_plot_max, max_mInc, max_mOrd, max_mQuanta, max_nInc, max_nInc2, max_nOrd, &
                     max_nOrd2, max_nQuanta, mTabDim, n, n_max2, n_plot_max, nTabDim, nTabDim2, nvTabDim
real(kind=wp) :: const1, rsum
integer(kind=iwp), allocatable :: mDec(:,:), mInc(:,:), mMat(:,:), nDec(:,:), nInc(:,:), nMat(:,:)
real(kind=wp), allocatable :: DipMat(:,:,:), FC2(:,:,:)

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
      rsum = Zero
      do j=1,nDimTot
        do i=1,nDimTot
          rsum = rsum+DipMat(i,j,k)*U1(i,m)*U2(j,n)
        end do
      end do
      FC2(m,n,k) = rsum
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

subroutine Intensity2(IntensityMat,TermMat,U1,U2,C1,W1,det1,r01,C2,W2,det2,r02,det0,m_max,n_max,max_dip,m_plot,n_plot,TranDip, &
                      TranDipGrad,Base,l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc,nDimTot,l_n_plot,l_m_plot)
!  Purpose:
!    Calculates the intensities of the different transitions between the two surfaces.
!
!  Input:
!    FC           : Real two dimensional array - Franck-Condon factors.
!    ipow         : Two dimensional integer array - terms of the polynomial.
!    var          : Real two dimensional array - coordinates to be used in the fit.
!    trfName      : Character array - transformation associated with each internal coordinate.
!    use_weight   : Logical
!    U1,U2        : Real two dimensional arrays - eigenvectors obtained from matrix element calculations.
!    W1,W2        : Real two dimensional arrays - eigenvectors scaled by the square root of the eigenvalues. Harmonic approximation.
!    C1,C2        : Real two dimensional arrays - inverses of W1 and W2.
!    det1,det2    : Real variables - determinants of C1 and C2.
!    r01,r02      : Real arrays - coordinates of the two oscillators.
!    m_max,n_max  : Integer variables - maximum quanta.
!    max_dip      : Integer variable - maximum order of transition dipole.
!    m_plot,
!    n_plot       : Integer array - transitions wanted in output.
!
!  Output:
!    IntensityMat : Real two dimensional array - intensities of the transitions.
!    TermMat      : Real two dimensional array - energies of transitions.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m_max, n_max, max_dip, l_n_plot, l_m_plot, m_plot(l_m_plot), n_plot(l_n_plot), l_IntensityMat_1, &
                                 l_IntensityMat_2, l_TermMat_1, l_TermMat_2, nOsc
real(kind=wp), intent(out) :: IntensityMat(0:l_IntensityMat_1,0:l_IntensityMat_2), det0
integer(kind=iwp), intent(inout) :: nDimTot
real(kind=wp), intent(in) :: TermMat(0:l_TermMat_1,0:l_TermMat_2), U1(nDimTot,nDimTot), U2(nDimTot,nDimTot), C1(nOsc,nOsc), &
                             W1(nOsc,nOsc), det1, r01(nOsc), C2(nOsc,nOsc), W2(nOsc,nOsc), det2, r02(nOsc), TranDip(3), &
                             TranDipGrad(3,nOsc), Base(nOsc,nOsc)
integer(kind=iwp) :: i, iOrd, j, jOrd, k, m, m_plot_max, max_mInc, max_mOrd, max_mQuanta, max_nInc, max_nInc2, max_nOrd, &
                     max_nOrd2, max_nQuanta, mTabDim, n, n_max2, n_plot_max, nTabDim, nTabDim2, nvTabDim
real(kind=wp) :: const1, rsum
integer(kind=iwp), allocatable :: mDec(:,:), mInc(:,:), mMat(:,:), nDec(:,:), nInc(:,:), nMat(:,:)
real(kind=wp), allocatable :: DipMat(:,:,:), FC2(:,:,:)

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
      rsum = Zero
      do j=1,nDimTot
        do i=1,nDimTot
          rsum = rsum+DipMat(i,j,k)*U1(i,m)*U2(j,n)
        end do
      end do
      FC2(m,n,k) = rsum
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

subroutine DipMatEl(Dij,W,L,U,FC00,nMat,nInc,nDec,D0,D1,D2,D3,D4,max_term,Base,m_ord,nOsc,max_mOrd,max_nOrd2)
!  Purpose:
!    Calculate matrix elements of the transition dipole.
!
!  Input:
!    D0       : Real variable - the zero order term of the transition dipole.
!    D1       : Real array - the first order term of the transition dipole.
!    D2       : Real two dimensional array - the second order term of the transition dipole.
!    D3       : Real three dimensional array - the third order term of the transition dipole.
!    D4       : Real four dimensional array - the fourth order term of the transition dipole.
!    W        : Real two dimensional array
!    L,U      : Real two dimensional array
!    FC00     : Real variable
!    nMat     : Two dimensional integer array.
!    max_term : Integer - maximum order of the transition dipole terms.
!
!  Output:
!    Dij      : Real two dimensional array - contains the matrix elements of the transition dipole.

use mula_global, only: ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2), max_term, m_ord, nOsc, max_mOrd, &
                                 max_nOrd2
real(kind=wp), intent(out) :: Dij(0:max_mOrd,0:max_mOrd)
real(kind=wp), intent(in) :: W(nOsc,nOsc), L(0:max_mOrd,0:max_mOrd), U(0:max_nOrd2,0:max_nord2), FC00, D0, D1(nOsc), &
                             D2(nOsc,nOsc), D3(nOsc,nOsc,nOsc), D4(nOsc,nOsc,nOsc,nOsc), Base(nOsc,nOsc)
integer(kind=iwp) :: max_nOrd, mPlus, nOscOld, nPlus
real(kind=wp), allocatable :: A(:,:), Temp(:,:), Wtemp(:,:)

! Initialize.
max_nOrd = max_mOrd
mPlus = max_mOrd+1
nPlus = max_nOrd+1
nOscOld = nOsc
call mma_allocate(A,[0,max_mOrd],[0,max_nOrd],label='A')
call mma_allocate(Wtemp,nOscOld,nOsc,label='Wtemp')
call DGEMM_('N','N',nOscOld,nOsc,nOsc,One,Base,nOscOld,W,nOsc,Zero,Wtemp,nOscOld)
A(:,:) = Zero
call PotEnergy(A,nMat,nInc,nDec,D0,D1,D2,D3,D4,max_term,WTemp,m_ord,nosc,nOscOld)

call mma_deallocate(Wtemp)
call mma_allocate(Temp,[0,max_mOrd],[0,max_nOrd],label='Temp')
call DGEMM_('N','T',mPlus,mPlus,nPlus,One,A,mPlus,U,mPlus,Zero,Temp,mPlus)
call DGEMM_('N','N',mPlus,nPlus,mPlus,FC00,L,mPlus,Temp(:,max_nOrd2+1:),mPlus,Zero,Dij,mPlus)
call mma_deallocate(Temp)
call mma_deallocate(A)

end subroutine DipMatEl

subroutine SetUpDipMat(DipMat,max_term,ipow,var,dip,trfName,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_mInc, &
                       max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,det0,r0,r1,r2,base,nOsc,nDimTot,nPolyTerm,ndata,nvar, &
                       MaxNumAt)
!  Purpose:
!    Performs a least squares fit of the transition dipole at the two
!    centers and at the inPolyTermediate oscillator. Calculates the matrix
!    elements of the transition dipole at these centers.
!
!  Input:
!    ipow       : Two dimensional integer array - terms of the polynomial.
!    var        : Real two dimensional array - coordinates to be used in the fit.
!    dip        : Real array - values of dipole at the coordinates contained in var.
!    trfName    : Character array - transformation associated with each internal coordinate.
!    max_term   : Integer - maximum order of the transition dipole terms.
!    W1,W2      : Real two dimensional arrays - eigenvectors scaled by the square root of the eigenvalues.
!    C1,C2      : Real two dimensional arrays - inverses of W1 and W2.
!    det1,det2  : Real variables - determinants of C1 and C2.
!    r01,r02    : Real arrays - coordinates of the two oscillators.
!    max_mOrd,
!    max_nOrd,
!    max_nOrd2
!    max_mInc,
!    max_nInc,
!    max_nInc2  : Integer variables
!    mMat,nMat,
!    mInc,nInc,
!    mDec,nDec  : Two dimensional integer arrays
!
!  Output:
!    DipMat     : Real two dimensional array - contains the matrix elements of the transition dipole.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: max_term, nPolyTerm, nvar, ipow(nPolyTerm,nvar), max_mOrd, max_nOrd, max_nOrd2, max_mInc, &
                                 max_nInc, max_nInc2, mMat(0:mdim1,mdim2), nMat(0:ndim1,ndim2), mInc(0:mdim1,mdim2), &
                                 nInc(0:ndim1,ndim2), mDec(0:mdim1,mdim2), nDec(0:ndim1,ndim2), nOsc, nDimtot, ndata, MaxNumAt
real(kind=wp), intent(out) :: DipMat(0:nDimTot,0:nDimTot), det0, r0(nOsc), r1(nOsc), r2(nOsc)
real(kind=wp), intent(in) :: var(ndata,nvar), dip(ndata), C1(nOsc,nOsc), W1(nOsc,nOsc), det1, r01(nOsc), C2(nOsc,nOsc), &
                             W2(nOsc,nOsc), det2, r02(nOsc), Base(nOsc,nOsc)
character(len=80), intent(in) :: trfName(MaxNumAt)
integer(kind=iwp) :: iOrd, jOrd, l_C1, l_C2, l_r0, l_r1, l_r2
real(kind=wp) :: D0, FC00, max_err, stand_dev
logical(kind=iwp) :: find_minimum, use_weight
real(kind=wp), allocatable :: alpha1(:,:), alpha2(:,:), beta(:,:), C(:,:), coef(:), D1(:), D2(:,:), D3(:,:,:), D4(:,:,:,:), &
                              Dij(:,:), DijTrans(:,:), L(:,:), r0vec(:), Sij(:,:), U(:,:), W(:,:)

! Initialize.
find_minimum = .false.
use_weight = .false.

call mma_allocate(Dij,[0,max_mOrd],[0,max_mOrd],label='Dij')
call mma_allocate(DijTrans,[0,max_mOrd],[0,max_mOrd],label='DijTrans')
call mma_allocate(C,nOsc,nOsc,label='C')
call mma_allocate(W,nOsc,nOsc,label='W')
call mma_allocate(L,[0,max_mOrd],[0,max_mOrd],label='L')
call mma_allocate(U,[0,max_nOrd2],[0,max_nOrd2],label='U')
call mma_allocate(Sij,[0,max_mOrd],[0,max_mOrd],label='Sij')
call mma_allocate(r0vec,nOsc,label='r0vec')
call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')
call mma_allocate(beta,nOsc,nOsc,label='beta')
call mma_allocate(D1,nOsc,label='D1')
call mma_allocate(D2,nOsc,nOsc,label='D2')
call mma_allocate(D3,nOsc,nOsc,nOsc,label='D3')
call mma_allocate(D4,nOsc,nOsc,nOsc,nOsc,label='D4')

! Calculate terms of type
!
!                    dM  |
!                    --  | < i | Q  | j >,
!                    dQ  |        k
!                      k  Q
!                          0
! where M is the (transition) dipole moment, |i> and |j> are harmonic
! oscillator states, Q_0 is the equilibrium geometry and Q_k is the
! k:th normal coordinate.

l_C1 = nOsc
call Calc_r00(C1,C1,C,W,alpha1,alpha2,r0vec,r01,r01,det0,det1,det1,FC00,l_C1)
call FCval(C1,W1,det1,r01,C1,W1,det1,r01,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det1,L,U,FC00,alpha1,alpha2,beta,l_C1)
call mma_allocate(coef,nPolyTerm,label='coef')
l_r1 = nOsc
call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,coef,r1,l_r1,D0,D1,D2,D3,D4,trfName,stand_dev,max_err,find_minimum,max_term, &
            use_weight,nOsc,nOsc,nOsc)
call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
DipMat(0:max_mOrd,0:max_mOrd) = Dij

l_C2 = nOsc
call Calc_r00(C2,C2,C,W,alpha1,alpha2,r0vec,r02,r02,det0,det2,det2,FC00,l_C2)
call FCval(C2,W2,det2,r02,C2,W2,det2,r02,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det2,L,U,FC00,alpha1,alpha2,beta,l_C2)
l_r2 = nOsc
call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,coef,r2,l_r2,D0,D1,D2,D3,D4,trfName,stand_dev,max_err,find_minimum,max_term, &
            use_weight,nOsc,nOsc,nOsc)
call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
DipMat(max_mOrd+1:2*max_mOrd+1,max_mOrd+1:2*max_mOrd+1) = Dij

l_C1 = nOsc
call Calc_r00(C1,C2,C,W,alpha1,alpha2,r0vec,r01,r02,det0,det1,det2,FC00,l_C1)
call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det0,L,U,FC00,alpha1,alpha2,beta,l_C1)
l_r0 = nOsc
call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,coef,r0,l_r0,D0,D1,D2,D3,D4,trfName,stand_dev,max_err,find_minimum,max_term, &
            use_weight,nOsc,nOsc,nOsc)
call mma_deallocate(coef)
call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
DipMat(0:max_mOrd,max_mOrd+1:2*max_mOrd+1) = Dij
do iOrd=0,max_mOrd
  do jOrd=0,max_mOrd
    DijTrans(jOrd,iOrd) = Dij(iOrd,jOrd)
  end do
end do
DipMat(max_mOrd+1:2*max_mOrd+1,0:max_mOrd) = DijTrans

call mma_deallocate(Dij)
call mma_deallocate(DijTrans)
call mma_deallocate(C)
call mma_deallocate(W)
call mma_deallocate(L)
call mma_deallocate(U)
call mma_deallocate(Sij)
call mma_deallocate(r0vec)
call mma_deallocate(alpha1)
call mma_deallocate(alpha2)
call mma_deallocate(beta)
call mma_deallocate(D1)
call mma_deallocate(D2)
call mma_deallocate(D3)
call mma_deallocate(D4)

end subroutine SetUpDipMat

!end module IntensityMod
