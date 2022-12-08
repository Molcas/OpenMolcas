!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Calc_r00(C1,C2,C,W,alpha1,alpha2,r00,r01,r02,det0,det1,det2,FC00,nOsc)
!  Purpose:
!    Calculate geometry of the intermediate oscillator.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOsc
real(kind=wp), intent(in) :: C1(nOsc,nOsc), C2(nOsc,nOsc), r01(nOsc), r02(nOsc), det1, det2
real(kind=wp), intent(out) :: C(nOsc,nOsc), W(nOsc,nOsc), alpha1(nOsc,nOsc), alpha2(nOsc,nOsc), r00(nOsc), det0, FC00
integer(kind=iwp) :: my1, my2
real(kind=wp) :: det, FC00_exp
real(kind=wp), allocatable :: alpha(:,:), beta(:,:), r_temp(:), r_temp1(:), r_temp2(:), temp(:,:), temp1(:,:)
real(kind=wp), external :: Ddot_

! Initialize.
my1 = nOsc
my2 = nOsc

call mma_allocate(temp,nOsc,nOsc,label='temp')

! Calculate alpha1, alpha2 and alpha.
call mma_allocate(alpha,nOsc,nOsc,label='alpha')
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C1,nOsc,C1,nOsc,Zero,alpha1,nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C2,nOsc,C2,nOsc,Zero,alpha2,nOsc)

temp(:,:) = alpha1+alpha2

alpha(:,:) = Half*temp

! Calculate C using a Cholesky factorization of 2*alpha.
call Cholesky(temp,C,nOsc)

! Calculate W.
call unitmat(W,nOsc)
temp(:,:) = C
call Dool_MULA(temp,nOsc,nOsc,W,my1,my2,det0)
det0 = abs(det0)

! Calculate r00.
call mma_allocate(r_temp1,nOsc,label='r_temp1')
call mma_allocate(r_temp2,nOsc,label='r_temp2')
call mma_allocate(r_temp,nOsc,label='r_temp')
call DGEMM_('N','N',nOsc,1,nOsc,One,alpha1,nOsc,r01,nOsc,Zero,r_temp1,nOsc)
call DGEMM_('N','N',nOsc,1,nOsc,One,alpha2,nOsc,r02,nOsc,Zero,r_temp2,nOsc)
r_temp(:) = r_temp1+r_temp2
temp(:,:) = Two*alpha

call Dool_MULA(temp,nOsc,nOsc,r_temp,nOsc,1,det)
call mma_allocate(beta,nOsc,nOsc,label='beta')
r00(:) = r_temp

! Calculate beta.
call mma_allocate(temp1,nOsc,nOsc,label='temp1')
temp1(:,:) = alpha1

temp(:,:) = Two*alpha

call Dool_MULA(temp,nOsc,nOsc,temp1,nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,alpha2,nOsc,temp1,nOsc,Zero,beta,nOsc)

! Calculate FC00.
r_temp1(:) = r01-r02

call DGEMM_('N','N',nOsc,1,nOsc,One,beta,nOsc,r_temp1,nOsc,Zero,r_temp2,nOsc)
FC00_exp = Ddot_(nOsc,r_temp1,1,r_temp2,1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)

call mma_deallocate(r_temp1)
call mma_deallocate(r_temp2)
call mma_deallocate(r_temp)
call mma_deallocate(alpha)
call mma_deallocate(beta)
call mma_deallocate(temp)
call mma_deallocate(temp1)

end subroutine Calc_r00
