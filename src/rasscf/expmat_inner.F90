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
! Jie J. Bao, on Feb. 03, 2022, created this file.               *
! ****************************************************************

subroutine ExpMat_Inner(R,X,nLen)
! Purpose: calculate R=exp(X)
! Explanation:
! The subroutine uses the algorithm in Section 3.1.5 on Page 83 of
! 'Molecular Electronic-Structure Theory' by T. Helgaker, P.
! Jorgensen and J. Olsen.
! Jie Bao acknowledges Dr. Chen Zhou from Xiamen University, China, for
! showing the resource of this algorithm and his code for reference.

use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer nLen, nLen2
real*8 x(nLen**2)
real*8 R(nLen**2)
real*8 tau2(nLen), cospart(nLen**2), sinpart(nLen**2), scr(nLen**2), X2(nLen**2), tau(nLen)
integer nScrDiag, INFO, I
real*8, dimension(:), allocatable :: ScrDiag
real*8 Coeff

nLen2 = nLen**2

!Step 1 calculate X2
call DGEMM_('n','n',nLen,nLen,nLen,1.0d0,X,nLen,X,nLen,0.0d0,X2,nLen)

!Step 2 diagonalize X2
call GetDiagScr(nScrDiag,X2,Scr,nLen)
call mma_allocate(ScrDiag,nScrDiag)
call DSYEV_('V','U',nLen,X2,nLen,tau2,ScrDiag,nScrDiag,INFO)
call mma_deallocate(ScrDiag)

do I=1,nLen
  tau(I) = dsqrt(dabs(tau2(I)))
end do

!Step 3 build cos part of R matrix
call DCopy_(nLen2,X2,1,CosPart,1)
do I=1,nLen
  call DScal_(nLen,cos(tau(I)),CosPart((I-1)*nLen+1),1)
end do

call DGEMM_('n','t',nLen,nLen,nLen,1.0d0,CosPart,nLen,X2,nLen,0.0d0,Scr,nLen)
! R = W * cos(tau) * W^T
call DCopy_(nLen2,Scr,1,R,1)
!Step 4 build sin part of R matrix
call DCopy_(nLen2,X2,1,SinPart,1)
do I=1,nLen
  if (tau(I) < 1.0d-8) then
    Coeff = 1.0d0
  else
    Coeff = sin(tau(I))/tau(I)
  end if
  call DScal_(nLen,Coeff,SinPart((I-1)*nLen+1),1)
end do

call DGEMM_('n','t',nLen,nLen,nLen,1.0d0,SinPart,nLen,X2,nLen,0.0d0,Scr,nLen)
! R  = R + W * tau^(-1) * sin(tau) * W^T * X
call DGEMM_('n','n',nLen,nLen,nLen,1.0d0,Scr,nLen,X,nLen,1.0d0,R,nLen)

return

end subroutine ExpMat_Inner
