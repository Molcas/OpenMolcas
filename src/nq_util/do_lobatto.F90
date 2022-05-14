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

! This subroutine should be in a module, to avoid explicit interfaces
#ifdef _IN_MODULE_

subroutine Do_Lobatto(L_Eff,nPoints,R)
!***********************************************************************
!                                                                      *
!     Computes data useful for the angular quadrature.                 *
!                                                                      *
!***********************************************************************

use nq_Grid, only: Pax
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L_Eff
integer(kind=iwp), intent(out) :: nPoints
real(kind=wp), allocatable, intent(out) :: R(:,:)
integer(kind=iwp) :: iOff, iOffT, iPhi, iTheta, mTheta, nLabatto, nPhi, nTheta
real(kind=wp) :: Cos_Phi, Cos_Theta, Sin_Phi, Sin_Theta, w_Phi, w_Theta, x, y, z
real(kind=wp), allocatable :: Labatto(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate angular grid a la Lobatto

nPoints = 0
nTheta = (L_Eff+3)/2
do iTheta=1,nTheta
  nPhi = L_Eff
  if ((iTheta == 1) .or. (iTheta == nTheta)) nPhi = 1
  if ((iTheta == nTheta/2+1) .and. (nTheta/2*2-nTheta == -1) .and. (nTheta > 3)) nPhi = L_Eff+4

  nPoints = nPoints+nPhi

end do

call mma_allocate(R,4,nPoints,Label='R')

nTheta = (L_Eff+3)/2
nLabatto = 3*(nTheta+2)*(nTheta+3)/2
call mma_allocate(Labatto,nLabatto,Label='Labatto')
call Lobatto(ntheta,Labatto)

mTheta = nTheta-1
iOffT = 1+3*mTheta*(mTheta+1)/2
iOff = 1
do iTheta=1,nTheta

  Cos_Theta = Labatto(iOffT)
  Sin_Theta = sqrt(One-Cos_Theta**2)
  w_Theta = Labatto(iOffT+1)
  iOffT = iOffT+3

  nPhi = L_Eff
  if ((iTheta == 1) .or. (iTheta == nTheta)) nPhi = 1
  if ((iTheta == nTheta/2+1) .and. (nTheta/2*2-nTheta == -1) .and. (nTheta > 3)) nPhi = L_Eff+4

  do iPhi=1,nPhi
    call Phi_point(iPhi,nPhi,Cos_Phi,Sin_Phi,w_Phi)

    x = Sin_Theta*Cos_Phi
    y = Sin_Theta*Sin_Phi
    z = Cos_Theta
    R(1,iOff) = Pax(1,1)*x+Pax(1,2)*y+Pax(1,3)*z
    R(2,iOff) = Pax(2,1)*x+Pax(2,2)*y+Pax(2,3)*z
    R(3,iOff) = Pax(3,1)*x+Pax(3,2)*y+Pax(3,3)*z
    R(4,iOff) = w_Theta*w_Phi
    iOff = iOff+1

  end do ! iPhi

end do   ! iTheta
call mma_deallocate(Labatto)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Do_Lobatto

#endif
