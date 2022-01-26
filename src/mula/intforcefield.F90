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

!contains

subroutine IntForceField(IntensityMat,TermMat,FC00,C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,m_max,n_max,max_dip,Trandip,TranDipGrad, &
                         harmfreq1,x_anharm1,harmfreq2,x_anharm2,mMat,mInc,mDec,nMat,nInc,nDec,OscStr,max_mOrd,max_nOrd,nOsc)
!  Purpose:
!    Calculates the intensities of the different transitions between the two surfaces.
!
!  Input:
!    FC           : Real two dimensional array - Franck-Condon factors.
!    W1,W2        : Real two dimensional arrays - eigenvectors scaled by the square root of the eigenvalues. Harmonic approximation.
!    C1,C2        : Real two dimensional arrays - inverses of W1 and W2.
!    det1,det2    : Real variables - determinants of C1 and C2.
!    r01,r02      : Real arrays - coordinates of the two oscillators.
!    m_max,n_max  : Integer variables - maximum quanta.
!    max_dip      : Integer variable - maximum order of transition dipole.
!    MatEl        : Logical
!    m_plot,
!    n_plot       : Integer array - transitions wanted in output.
!    TermMat      : Real two dimensional array - energies of transitions.
!
!  Output:
!    IntensityMat : Real two dimensional array - intensities of the transitions.

use mula_global, only: mdim1, mdim2, ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, auTokJ, kBoltzmann
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: max_mOrd, max_nOrd
real(kind=wp), intent(out) :: IntensityMat(0:max_mOrd,0:max_nOrd), FC00
integer(kind=iwp), intent(in) :: m_max, n_max, max_dip, mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2), &
                                 nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2), nOsc
real(kind=wp), intent(in) :: TermMat(0:max_mOrd,0:max_nOrd), C1(nOsc,nOsc), W1(nOsc,nOsc), det1, r01(nOsc), C2(nOsc,nOsc), &
                             W2(nOsc,nOsc), det2, r02(nOsc), C(nOsc,nOsc), W(nOsc,nOsc), det0, TranDip(3), TranDipGrad(3,nOsc), &
                             harmfreq1(nOsc), x_anharm1(nOsc,nOsc), harmfreq2(nOsc), x_anharm2(nOsc,nOsc)
logical(kind=iwp), intent(in) :: OscStr
integer(kind=iwp) :: iOrd, jOrd, l_harm, nvTabDim
real(kind=wp) :: const1, const2, dE
integer(kind=iwp), allocatable :: level1(:), level2(:)
real(kind=wp), allocatable :: FC2(:,:,:), FreqDiffmat(:)
real(kind=wp), parameter :: Temperature = 10.0_wp

! Initialize.
call TabDim(m_max,nosc,nvTabDim)
max_mOrd = nvTabDim-1
call TabDim(n_max,nosc,nvTabDim)
max_nOrd = nvTabDim-1
call mma_allocate(FC2,[0,max_mOrd],[0,max_nOrd],[0,3],label='FC2')

call SetUpHarmDip(FC2,max_dip,m_max,n_max,mMat,mInc,mDec,nMat,nInc,nDec,C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,TranDip, &
                  TranDipGrad,FC00,max_mOrd,max_nOrd,nOsc)

! Calculate intensities with Boltzmann weighting of hotband intensity.
const1 = Two/Three
call mma_allocate(level1,nOsc,label='level1')
call mma_allocate(level2,nOsc,label='level2')

if (max_nOrd > max_mOrd) then
  call mma_allocate(FreqDiffmat,[0,max_mOrd],label='FreqDiffMat')

  level1(:) = mMat(0,:)
  do iOrd=0,max_mOrd
    level2(:) = mMat(iOrd,:)
    l_harm = nOsc
    call TransEnergy(Zero,x_anharm1,harmfreq1,level1,Zero,x_anharm1,harmfreq1,level2,FreqDiffMat(iOrd),l_harm)
  end do
  do jOrd=0,max_nOrd
    do iOrd=0,max_mOrd
      dE = FreqDiffMat(iOrd)*auTokJ*1.0e-39_wp
      const2 = const1*exp(-dE/(kBoltzmann*Temperature))
      IntensityMat(iOrd,jOrd) = const2*abs(TermMat(iOrd,jOrd))*(FC2(iOrd,jOrd,1)**2+FC2(iOrd,jOrd,2)**2+FC2(iOrd,jOrd,3)**2)
      if (.not. OscStr) then
        ! where does this number come from?
        IntensityMat(iOrd,jOrd) = 32.13002e9_wp*const2*(TermMat(iOrd,jOrd)**2)*IntensityMat(iOrd,jOrd)
      end if
    end do
  end do
else
  call mma_allocate(FreqDiffMat,[0,max_nOrd],label='FreqDiffMat')
  level1(:) = nMat(0,:)
  do iOrd=0,max_nOrd
    level2(:) = nMat(iOrd,:)
    l_harm = nOsc
    call TransEnergy(Zero,x_anharm2,harmfreq2,level1,Zero,x_anharm2,harmfreq2,level2,FreqDiffMat(iOrd),l_harm)
  end do
  do jOrd=0,max_nOrd
    dE = FreqDiffMat(jOrd)*auTokJ*1.0e-39_wp
    do iOrd=0,max_mOrd
      const2 = const1*exp(-dE/(kBoltzmann*Temperature))
      IntensityMat(iOrd,jOrd) = const2*abs(TermMat(iOrd,jOrd))*(FC2(iOrd,jOrd,1)**2+FC2(iOrd,jOrd,2)**2+FC2(iOrd,jOrd,3)**2)
      if (.not. OscStr) then
        ! where does this number come from?
        IntensityMat(iOrd,jOrd) = 32.13002e9_wp*const2*(TermMat(iOrd,jOrd)**2)*IntensityMat(iOrd,jOrd)
      end if
    end do
  end do
end if
call mma_deallocate(level1)
call mma_deallocate(level2)
call mma_deallocate(FC2)
call mma_deallocate(FreqDiffMat)

end subroutine IntForceField

!end module ForceFieldIntMod
