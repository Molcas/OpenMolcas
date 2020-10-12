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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
SUBROUTINE pgek()
use kriging_mod, only: x, y, nPoints, nInter_Save
Implicit None
#include "real.fh"
#include "stdalloc.fh"
Real*8, Allocatable:: Mean_univariate(:), Mean_bivariate(:)
Real*8, Allocatable:: Variance_univariate(:), Variance_bivariate(:)
Integer i, j, k, l
Real*8 tmp, dx, dy
! Mutual information array
!Real*8 MI(nInter_save)
Real*8, Allocatable::  MI(:)
! universal kernal density estimators.
!Real*8 px(nPoints,nInter_save)
!Real*8 py(nPoints)
!Real*8 pxy(nPoints,nInter_save)
Real*8, Allocatable:: px(:,:), py(:), pxy(:,:)
! the 2 x 2 covariance matrix (used to compute the transformation variable in the bivariate case)
Real*8 Sigma(2,2), Sigma_Inverse(2,2)
Real*8 d, d_h, h, u_j, K_u_j, Norm_Sigma

Write (6,*) 'PGEK: nPoints=',nPoints
Write (6,*) 'PGEK: nInter_Save=',nInter_Save
Call RecPrt('y',' ',y,1,nPoints)
Call RecPrt('x',' ',x,nInter_Save,nPoints)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    Allocate memory

Call mma_allocate(Mean_univariate,nInter_save+1,Label='Mean_univariate')
Call mma_allocate(Mean_bivariate,nInter_save,Label='Mean_bivariate')
Call mma_allocate(Variance_univariate,nInter_save+1,Label='Variance_univariate')
Call mma_allocate(Variance_bivariate,nInter_save,Label='Variance_bivariate')
Call mma_allocate(MI,nInter_Save,Label='MI')
Call mma_allocate(px,nPoints,nInter_Save,Label='px')
Call mma_allocate(py,nPoints,Label='py')
Call mma_allocate(pxy,nPoints,nInter_Save,Label='pxy')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Compute the mean for all variables - coordinate components and the energy
!  univariate --  x and y -- and bivariate -- xy
!  Compute the variance of the same

Do i = 1, nInter_save
  Mean_Univariate(i)    =Zero
  Mean_bivariate(i)     =Zero
  Variance_Univariate(i)=Zero
  Variance_bivariate(i) =Zero
  Do j = 1, nPoints
    Mean_Univariate(i)    =Mean_Univariate(i)    + x(i,j)
    Mean_bivariate(i)     =Mean_bivariate(i)     + x(i,j)*y(j)
    Variance_Univariate(i)=Variance_Univariate(i)+ (x(i,j))**2
    Variance_bivariate(i) =Variance_bivariate(i) + (x(i,j)*y(j))**2
  End Do
End Do
Mean_Univariate(nInter_save+1)    =Zero
Variance_Univariate(nInter_save+1)=Zero
Do j = 1, nPoints
  Mean_Univariate(nInter_save+1)    =Mean_Univariate(nInter_save+1)     + y(j)
  Variance_Univariate(nInter_save+1)=Variance_Univariate(nInter_save+1) + y(j)**2
End Do
Mean_Univariate(:)    = Mean_Univariate(:)    /DBLE(nPoints)
Mean_Bivariate(:)     = Mean_Bivariate(:)     /DBLE(nPoints)
Variance_Univariate(:)= Variance_Univariate(:)/DBLE(nPoints)
Variance_Bivariate(:) = Variance_Bivariate(:) /DBLE(nPoints)

Call RecPrt('Mean - x ',' ',Mean_Univariate(:),nInter_Save,1)
Call RecPrt('Mean - xy ',' ',Mean_Bivariate(:),nInter_Save,1)
Call RecPrt('Variance - x ',' ',Variance_Univariate(:),nInter_Save,1)
Call RecPrt('Variance - xy ',' ',Variance_Bivariate(:),nInter_Save,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the Gaussian bandwidth

d_h = Two
h = (Four/(d_h + Two))**(One/(Four+d_h)) * DBLE(nPoints)**(-One/(Four+d_h))
Write (6,*) 'd_h=',d_h
Write (6,*) 'h=',h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the universal kernel density estimators px, py, and pxy
! also known as Parzen window

! 1) the univariate cases

d=One
Write (6,*) 'd=',d


!  a) Compute px(i,l)

Do l = 1, nInter_save
   Do i = 1, nPoints

      Norm_Sigma=Abs(Variance_Univariate(l))
      tmp=Zero
      Do j = 1, nPoints
         u_j = (x(l,i)-x(l,j))**2 / (h**2 * Variance_univariate(l))
         K_u_j = One/((Two*Pi)**(d/Two) * h**d * sqrt(Norm_Sigma)) * Exp(-u_j/Two)
         tmp=tmp+K_u_j
      End Do

      px(i,l)=tmp/DBLE(nPoints)

   End Do
End Do
Call RecPrt('px',' ',px,nPoints,nInter_Save)

!  b) Compute py(i)

Do i = 1, nPoints

   Norm_Sigma=Abs(Variance_Univariate(l))
   tmp=Zero
   Do j = 1, nPoint
      u_j = (y(i)-y(j))**2 / (h**2 * Variance_univariate(l))
      K_u_j = One/((Two*Pi)**(d/Two) * h**d * sqrt(Norm_Sigma)) * Exp(-u_j/Two)
      tmp=tmp+K_u_j
   End Do

   py(i)=tmp/DBLE(nPoints)

End Do
Call RecPrt('py',' ',py,nPoints,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2) the bivariate case, pxy(i,l)

d=Two
Write (6,*) 'd=',d

Do l = 1, nInter_save
   Do i = 1, nPoints

!     Assemble the 2 x 2 covariance matrix
      Sigma(1,1)=Variance_univariate(l)  ! p(x_l^2)
      Sigma(1,2)=Variance_bivariate(l)   ! p(x_ly)
      Sigma(2,1)=Sigma(1,2)
      Sigma(2,2)=Variance_univariate(nInter_save+1) ! p(y^2)
      Call RecPrt('Sigma',' ',Sigma,2,2)
      Norm_Sigma=Sigma(1,1)*Sigma(2,2)-Sigma(1,2)*Sigma(2,1)  ! The determinant of the 2 x 2 sigma matrix
      Write (6,*) 'Norm_Sigma=',Norm_Sigma
      Sigma_Inverse(1,1)= Sigma(2,2)/Norm_Sigma
      Sigma_Inverse(1,2)=-Sigma(1,2)/Norm_Sigma
      Sigma_Inverse(2,1)=-Sigma(2,1)/Norm_Sigma
      Sigma_Inverse(2,2)=-Sigma(1,1)/Norm_Sigma
      Call RecPrt('Sigma_Inverse',' ',Sigma_Inverse,2,2)

      tmp=Zero
      Do j = 1, nPoints
         dx = x(l,i)-x(l,j)
         dy = y(i)-y(j)
         u_j = ( dx*dx * Sigma_Inverse(1,1)   &
                +dx*dy * Sigma_Inverse(1,2)   &
                +dy*dx * Sigma_Inverse(2,1)   &
                +dy*dy * Sigma_Inverse(2,2) ) &
               / h**2
         Write (6,*)' u_j=',u_j
         K_u_j = One/((Two*Pi)**(d/Two) * h**d * sqrt(Abs(Norm_Sigma))) * Exp(-u_j/Two)
         Write (6,*)' K_u_j=',K_u_j
         tmp=tmp+K_u_j
      End Do

      pxy(i,l)=tmp/DBLE(nPoints)

   End Do
End Do

Call RecPrt('pxy',' ',pxy,nPoints,nInter_Save)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assemble the Mutual Information, MI(l), for each variable l

Do l = 1, nInter_save
   tmp=0.0D0
   Do i = 1, nPoints
      tmp = tmp + pxy(i,l) * log10(pxy(i,l)/(px(i,l)*py(i)))
   End Do
   MI(l)=tmp
End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Write (*,*) 'Mutual Information'
Write (*,*) '=================='
Do i = 1, nInter_Save
   Write (6,'(i4,E10.2)') i, MI(i)
End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Deallocate memory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Call mma_deallocate(MI)
Call mma_deallocate(px)
Call mma_deallocate(py)
Call mma_deallocate(pxy)
Call mma_deallocate(Variance_univariate)
Call mma_deallocate(Variance_bivariate)
Call mma_deallocate(Mean_bivariate)
Call mma_deallocate(Mean_univariate)

END SUBROUTINE pgek
