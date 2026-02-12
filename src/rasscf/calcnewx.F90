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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 12, 2022, created this file.               *
! ****************************************************************

subroutine CalcNewX(X,H,G,nSPair,XScr,GScr,EigVal,ScrDiag,nScr)

! Commented lines are options under development and may be used in future
use CMS, only: CMSThres, PosHess, BigQaaGrad, nPosHess, LargestQaaGrad, NeedMoreStep

implicit none
integer nSPair, INFO, iPair, nScr
real*8 X(nSPair), G(nSPair), XScr(nSPair)
real*8 H(nSPair**2), ScrDiag(nScr), GScr(nSPair)
real*8 EigVal(nSPair**2)
real*8 MinGrad, ThreG, ThreH
real*8 ValGrad, AbsGrad, ValHess, AbsHess
! Thanks to Matthew R. Hermes for this algorithm

! Solve for x in hx=-g.
! 1. diagonalize h to h' with rotation matrix C
! 2. (C^T)h(C)(C^T)x=(C^T)g
! 3. define x'=(C^T)x; g'=(C^T)g, we have (h')(x')=(g')
! 4. (x')_p = (g')_p / (h')_pp
! 5. x=(C)(C^T)x=(C)(x')

!write(6,*) 'gradient before diag'
!call RecPrt(' ',' ',G,1,nSPair)
!write(6,*) 'hessian before diag'
!call RecPrt(' ',' ',EigVal,1,nSPair)
!*****Step 1
call DSYEV_('V','U',nSPair,H,nSPair,EigVal,ScrDiag,nScr,INFO)

!*****Step 3, g'=(C^T)g
call DGEMM_('n','n',1,nSPair,nSPair,1.0d0,G,1,H,nSPair,0.0d0,GScr,1)

!*****Step 4
ThreG = CMSThres*1.0d-2
ThreH = CMSThres*1.0d-4
MinGrad = CMSThres*1.0d5
!write(6,*) 'gradient'
!call RecPrt(' ','(1X,15(F9.6,1X))',GScr,1,nSPair)
!write(6,*) 'hessian'
!call RecPrt(' ','(1X,15(F9.6,1X))',EigVal,1,nSPair)

LargestQaaGrad = 0.0d0
nPosHess = 0
do iPair=1,nSPair
  ValGrad = GScr(iPair)
  AbsGrad = abs(ValGrad)
  ValHess = EigVal(iPair)
  AbsHess = abs(ValHess)

  if (ValHess > ThreH) nPosHess = nPosHess+1
  if (AbsGrad > LargestQaaGrad) LargestQaaGrad = AbsGrad

  if ((AbsGrad < ThreG) .and. (AbsHess < ThreH)) then
    !write(6,*) 'constant Qaa for pair',ipair
    XScr(iPair) = 0.0d0
  else
    XScr(iPair) = ValGrad/AbsHess
  end if

  if ((AbsGrad < ThreG) .and. (ValHess > ThreH)) then
    !write(6,*) 'local minimum for pair',ipair
    XScr(iPair) = MinGrad/abs(EigVal(iPair))
    if (XScr(iPair) > 1.0d-2) XScr(iPair) = 1.0d-2
  end if

end do

PosHess = .false.
BigQaaGrad = .false.
NeedMoreStep = .false.
if (nPosHess > 0) PosHess = .true.
if (LargestQaaGrad > ThreG) BigQaaGrad = .true.

if (PosHess .or. BigQaaGrad) NeedMoreStep = .true.

!write(6,*) 'steps taken'
!call RecPrt(' ','(1X,15(F9.6,1X))',XScr,1,nSPair)
!*****Step 5
call DGEMM_('n','t',1,nSPair,nSPair,1.0d0,XScr,1,H,nSPair,0.0d0,X,1)

return

end subroutine CalcNewX
