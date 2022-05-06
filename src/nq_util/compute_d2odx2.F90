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

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 ZA(nAtoms), O(3,3), EVal(3), dMdx(3,3), Px(3,3), dMdy(3,3), Py(3,3), d2Odx2(3,3)
logical Rot_Corr
! Local Arrays
real*8 d2Mdx2(3,3), Pxy(3,3), RHS(3,3), Scr1(3,3), Scr2(3,3), Scr3(3,3)

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Rot_Corr) then
  call FZero(d2Odx2,9)
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

call FZero(RHS,9)

! - O^T M^(xy) O

call DGEMM_('T','N',3,3,3,1.0d0,O,3,d2Mdx2,3,0.0d0,Scr1,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr1,3,O,3,0.0d0,Scr2,3)
call DaXpY_(9,-One,Scr2,1,RHS,1)

call FZero(Scr3,9)
Scr3(1,1) = Eval(1)
Scr3(2,2) = Eval(2)
Scr3(3,3) = Eval(3)

! + P^x Lambda P^y

call DGEMM_('N','N',3,3,3,1.0d0,Px,3,Scr3,3,0.0d0,Scr1,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr1,3,Py,3,0.0d0,Scr2,3)
call DaXpY_(9,One,Scr2,1,RHS,1)

! + P^y Lambda P^x

call DGEMM_('N','N',3,3,3,1.0d0,Py,3,Scr3,3,0.0d0,Scr1,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr1,3,Px,3,0.0d0,Scr2,3)
call DaXpY_(9,One,Scr2,1,RHS,1)

! + P^x O^T M^y O

call DGEMM_('N','T',3,3,3,1.0d0,Px,3,O,3,0.0d0,Scr1,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr1,3,dMdy,3,0.0d0,Scr2,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr2,3,O,3,0.0d0,Scr1,3)
call DaXpY_(9,One,Scr1,1,RHS,1)

! + P^y O^T M^x O

call DGEMM_('N','T',3,3,3,1.0d0,Py,3,O,3,0.0d0,Scr1,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr1,3,dMdx,3,0.0d0,Scr2,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr2,3,O,3,0.0d0,Scr1,3)
call DaXpY_(9,One,Scr1,1,RHS,1)

! - O^T M^x O P^y

call DGEMM_('T','N',3,3,3,1.0d0,O,3,dMdx,3,0.0d0,Scr1,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr1,3,O,3,0.0d0,Scr2,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr2,3,Py,3,0.0d0,Scr1,3)
call DaXpY_(9,-One,Scr1,1,RHS,1)

! - O^T M^y O P^x

call DGEMM_('T','N',3,3,3,1.0d0,O,3,dMdy,3,0.0d0,Scr1,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr1,3,O,3,0.0d0,Scr2,3)
call DGEMM_('N','N',3,3,3,1.0d0,Scr2,3,Px,3,0.0d0,Scr1,3)
call DaXpY_(9,-One,Scr1,1,RHS,1)
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

call DGEMM_('N','N',3,3,3,1.0d0,O,3,Pxy,3,0.0d0,d2Odx2,3)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Compute_d2Odx2
