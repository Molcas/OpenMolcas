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

subroutine Compute_d2Odx2(ZA,nAtoms,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,Px,jAtom,jCar,dMdy,Py,d2Odx2)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, iAtom, iCar, jAtom, jCar
real(kind=wp), intent(in) :: ZA(nAtoms), O(3,3), EVal(3), dTdRAi, dMdx(3,3), Px(3,3), dMdy(3,3), Py(3,3)
logical(kind=iwp), intent(in) :: Rot_Corr
real(kind=wp), intent(out) :: d2Odx2(3,3)
real(kind=wp) :: Alpha1, Alpha2, Beta1, Beta2, c12, c13, c23, d2Mdx2(3,3), Gamma1, Gamma2, Pxy(3,3), RHS(3,3), Scr1(3,3), &
                 Scr2(3,3), Scr3(3,3)

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Rot_Corr) then
  d2Odx2(:,:) = Zero
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!
!     Compute d2M/dxdy
!
call Compute_d2Mdx2(ZA,nAtoms,iAtom,iCar,dTdRAi,jAtom,jCar,d2Mdx2)

!                                                                      *
!***********************************************************************
!                                                                      *
!     Form diagonal elements of Pxy directly from Px and Py
!
Gamma1 = Px(1,2)
Beta1 = Px(3,1)
Alpha1 = Px(2,3)
Gamma2 = Py(1,2)
Beta2 = Py(3,1)
Alpha2 = Py(2,3)
Pxy(1,1) = -Gamma1*Gamma2-Beta1*Beta2
Pxy(2,2) = -Gamma1*Gamma2-Alpha1*Alpha2
Pxy(3,3) = -Beta1*Beta2-Alpha1*Alpha2
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute additional constraints for off-diagonals
!
c12 = Beta1*Alpha2+Beta2*Alpha1   ! Pxy(2,1)+Pxy(1,2)
c13 = Gamma1*Alpha2+Gamma2*Alpha1 ! Pxy(3,1)+Pxy(1,3)
c23 = Beta1*Gamma2+Beta2*Gamma1   ! Pxy(2,3)+Pxy(3,2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Assemble the right hand side of eq. 26 except for Lambda^(xy).
! This will generate all off-diagonal elements of the RHS!

! - O^T M^(xy) O

call DGEMM_('T','N',3,3,3,One,O,3,d2Mdx2,3,Zero,Scr1,3)
call DGEMM_('N','N',3,3,3,One,Scr1,3,O,3,Zero,Scr2,3)
RHS(:,:) = -Scr2

Scr3(:,:) = Zero
Scr3(1,1) = Eval(1)
Scr3(2,2) = Eval(2)
Scr3(3,3) = Eval(3)

! + P^x Lambda P^y

call DGEMM_('N','N',3,3,3,One,Px,3,Scr3,3,Zero,Scr1,3)
call DGEMM_('N','N',3,3,3,One,Scr1,3,Py,3,Zero,Scr2,3)
RHS(:,:) = RHS+Scr2

! + P^y Lambda P^x

call DGEMM_('N','N',3,3,3,One,Py,3,Scr3,3,Zero,Scr1,3)
call DGEMM_('N','N',3,3,3,One,Scr1,3,Px,3,Zero,Scr2,3)
RHS(:,:) = RHS+Scr2

! + P^x O^T M^y O

call DGEMM_('N','T',3,3,3,One,Px,3,O,3,Zero,Scr1,3)
call DGEMM_('N','N',3,3,3,One,Scr1,3,dMdy,3,Zero,Scr2,3)
call DGEMM_('N','N',3,3,3,One,Scr2,3,O,3,Zero,Scr1,3)
RHS(:,:) = RHS+Scr1

! + P^y O^T M^x O

call DGEMM_('N','T',3,3,3,One,Py,3,O,3,Zero,Scr1,3)
call DGEMM_('N','N',3,3,3,One,Scr1,3,dMdx,3,Zero,Scr2,3)
call DGEMM_('N','N',3,3,3,One,Scr2,3,O,3,Zero,Scr1,3)
RHS(:,:) = RHS+Scr1

! - O^T M^x O P^y

call DGEMM_('T','N',3,3,3,One,O,3,dMdx,3,Zero,Scr1,3)
call DGEMM_('N','N',3,3,3,One,Scr1,3,O,3,Zero,Scr2,3)
call DGEMM_('N','N',3,3,3,One,Scr2,3,Py,3,Zero,Scr1,3)
RHS(:,:) = RHS-Scr1

! - O^T M^y O P^x

call DGEMM_('T','N',3,3,3,One,O,3,dMdy,3,Zero,Scr1,3)
call DGEMM_('N','N',3,3,3,One,Scr1,3,O,3,Zero,Scr2,3)
call DGEMM_('N','N',3,3,3,One,Scr2,3,Px,3,Zero,Scr1,3)
RHS(:,:) = RHS-Scr1
#ifdef _DEBUGPRINT_
call RecPrt('RHS',' ',RHS,3,3)
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the off-diagonal elements.

! We will need some more elaborate code if the denominator is
! degenerate! Will be developed later...

Pxy(2,1) = (RHS(1,2)-EVal(1)*c12)/(EVal(2)-EVal(1))
Pxy(1,2) = c12-Pxy(2,1)
Pxy(3,1) = (RHS(1,3)-EVal(1)*c13)/(EVal(3)-EVal(1))
Pxy(1,3) = c13-Pxy(3,1)
Pxy(3,2) = (RHS(2,3)-EVal(2)*c23)/(EVal(3)-EVal(2))
Pxy(2,3) = c23-Pxy(3,2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Finally for O^(xy) from O P^(xyz)

call DGEMM_('N','N',3,3,3,One,O,3,Pxy,3,Zero,d2Odx2,3)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Compute_d2Odx2
