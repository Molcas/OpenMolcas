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
#define _DEBUGPRINT_
SUBROUTINE pgek()
use kriging_mod, only: x, y, nPoints, nInter
Implicit None
#include "real.fh"
#include "stdalloc.fh"
Real*8, Allocatable:: Mean_univariate(:)
Real*8, Allocatable:: Variance_univariate(:), Variance_bivariate(:)
Integer i, j,  l
Real*8 tmp, dx, dy, Fact
! Mutual information array
!Real*8 MI(nInter)
Real*8, Allocatable::  MI(:)
! universal kernal density estimators.
!Real*8 px(nPoints,nInter)
!Real*8 py(nPoints)
!Real*8 pxy(nPoints,nInter)
Real*8, Allocatable:: px(:,:), py(:), pxy(:,:)
! the 2 x 2 covariance matrix (used to compute the transformation variable in the bivariate case)
Real*8 Sigma2(2,2), Sigma2_Inverse(2,2)
Real*8 d, d_h, h, u_j, K_u_j, Det_Sigma2, Sigma_2, N_norm

#ifdef _DEBUGPRINT_
Write (6,*) 'PGEK: nPoints=',nPoints
Write (6,*) 'PGEK: nInter=',nInter
Call RecPrt('y',' ',y,1,nPoints)
Call RecPrt('x',' ',x,nInter,nPoints)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    Allocate memory

Call mma_allocate(Mean_univariate,nInter+1,Label='Mean_univariate')
Call mma_allocate(Variance_univariate,nInter+1,Label='Variance_univariate')
Call mma_allocate(Variance_bivariate,nInter,Label='Variance_bivariate')
Call mma_allocate(MI,nInter,Label='MI')
Call mma_allocate(px,nPoints,nInter,Label='px')
Call mma_allocate(py,nPoints,Label='py')
Call mma_allocate(pxy,nPoints,nInter,Label='pxy')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Compute the mean for all variables - coordinate components and the energy
!  univariate --  x and y
!  Compute the variance of the univariate and the bivariate

! univariate means
Do i = 1, nInter
  Mean_Univariate(i)    =Zero
  Do j = 1, nPoints
    Mean_Univariate(i)    =Mean_Univariate(i)    + x(i,j)
  End Do
End Do
Mean_Univariate(nInter+1)    =Zero
Do j = 1, nPoints
  Mean_Univariate(nInter+1)    =Mean_Univariate(nInter+1)     + y(j)
End Do
Mean_Univariate(:)    = Mean_Univariate(:)    / DBLE(nPoints)

! uni- and bivariate variances
Do i = 1, nInter
  Variance_Univariate(i)=Zero
  Variance_bivariate(i) =Zero
  Do j = 1, nPoints
    Variance_Univariate(i)=Variance_Univariate(i)+ (x(i,j) - Mean_univariate(i))**2
    Variance_bivariate(i) =Variance_bivariate(i) + (x(i,j) - Mean_univariate(i))  &
                                                 * (y(j)   - Mean_univariate(nInter+1))
  End Do
End Do
Variance_Univariate(nInter+1)=Zero
Do j = 1, nPoints
  Variance_Univariate(nInter+1)=Variance_Univariate(nInter+1) + (y(j)-Mean_univariate(nInter+1))**2
End Do
Variance_Univariate(:)= Variance_Univariate(:)/DBLE(nPoints)
Variance_Bivariate(:) = Variance_Bivariate(:) /DBLE(nPoints)

#ifdef _DEBUGPRINT_
Call RecPrt('Mean - x,y ',' ',Mean_Univariate(:),nInter+1,1)
Call RecPrt('Variance - x^2,y^2 ',' ',Variance_Univariate(:),nInter+1,1)
Call RecPrt('Variance - xy ',' ',Variance_Bivariate(:),nInter,1)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the Gaussian bandwidth

d_h = Two
h = (Four/(d_h + Two))**(One/(Four+d_h)) * DBLE(nPoints)**(-One/(Four+d_h))
#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'd_h=',d_h
Write (6,*) 'h=',h
Write (6,*)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the universal kernel density estimators px, py, and pxy
! also known as Parzen window technique using as a kernel function a normal distribution with
! an exponent of alpha=1/(2 h^2 Sigma^2), in which h is the Gaussian bandwith and Sigma^2 is the
! variance of the random variable under consideration -- Sigma is the standard deviation.

! 1) the univariate cases

d=One
Fact = One /  ( (Two*Pi)**(d/Two) * h**d )  ! Part of the normalization factor
#ifdef _DEBUGPRINT_
Write (6,*) 'd=',d
Write (6,*) '1/((2pi)^(d/2) h^d)=',Fact

#endif


!  a) Compute px(i,l)

Do l = 1, nInter
   Sigma_2=Variance_Univariate(l)
   Write (6,*) 'l, Sigma_2=',l, Sigma_2
!
!  Take special care of coordinates which might not change due to symmetry
!  or for some other reason. Then the kernel function effectively is a Dirac delta function.
!
   If (Sigma_2<1.0D-10) Then
      px(:,l)=One/DBLE(nPoints)
   Else
      N_norm= (Fact / SQRT(Sigma_2))
      Do i = 1, nPoints
         Write (6,*) 'i=',i

         tmp=Zero
         Do j = 1, nPoints
            u_j = (x(l,i)-x(l,j))**2 / (h**2 * Sigma_2)
            K_u_j = N_norm * EXP(-u_j/Two)
            Write (6,*) 'u_j, K_u_j=',u_j, K_u_j
            tmp=tmp+K_u_j
         End Do

         px(i,l)=tmp/DBLE(nPoints)

      End Do
   End If
End Do

!  b) Compute py(i)

Sigma_2=Variance_Univariate(nPoints+1)
N_norm= (Fact / SQRT(Sigma_2))
Write (6,*) 'Sigma_2=',Sigma_2
Do i = 1, nPoints
   Write (6,*) 'i=',i

   tmp=Zero
   Do j = 1, nPoints
      u_j = (y(i)-y(j))**2 / (h**2 * Sigma_2)
      K_u_j = N_norm * EXP(-u_j/Two)
      Write (6,*) 'u_j, K_u_j=',u_j, K_u_j
      tmp=tmp+K_u_j
   End Do

   py(i)=tmp/DBLE(nPoints)

End Do
#ifdef _DEBUGPRINT_
Call RecPrt('px',' ',px,nPoints,nInter)
Call RecPrt('py',' ',py,nPoints,1)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2) the bivariate case, pxy(i,l)

d=Two
Fact = One / ( (Two*Pi)**(d/Two) * h**d )
Write (6,*) 'd=',d
Write (6,*) 'Fact=',Fact

Do l = 1, nInter
   Sigma_2=Variance_univariate(l)
   If (Sigma_2<1.0D-10) Then
      pxy(:,l)=One/DBLE(nPoints)  ! Not sure about this part yet.
   Else
!     Assemble the 2 x 2 covariance matrix
      Sigma2(1,1)=Variance_univariate(l)  ! p(x_l^2)
      Sigma2(1,2)=Variance_bivariate(l)   ! p(x_ly)
      Sigma2(2,1)=Sigma2(1,2)
      Sigma2(2,2)=Variance_univariate(nInter+1) ! p(y^2)
      Call RecPrt('Sigma2',' ',Sigma2,2,2)
      Det_Sigma2=Sigma2(1,1)*Sigma2(2,2)-Sigma2(1,2)*Sigma2(2,1)  ! The determinant of the 2 x 2 sigma matrix
      Write (6,*) 'Det(Sigma2)=',Det_Sigma2
      Sigma2_Inverse(1,1)= Sigma2(2,2)/Det_Sigma2
      Sigma2_Inverse(1,2)=-Sigma2(1,2)/Det_Sigma2
      Sigma2_Inverse(2,1)=-Sigma2(2,1)/Det_Sigma2
      Sigma2_Inverse(2,2)= Sigma2(1,1)/Det_Sigma2
      Call RecPrt('Sigma2_Inverse',' ',Sigma2_Inverse,2,2)

      N_norm = (Fact / SQRT(Det_Sigma2))
      Do i = 1, nPoints


         tmp=Zero
         Do j = 1, nPoints
            dx = x(l,i)-x(l,j)
            dy = y(i)-y(j)
            u_j = ( dx*dx * Sigma2_Inverse(1,1)   &
                   +dx*dy * Sigma2_Inverse(1,2)   &
                   +dy*dx * Sigma2_Inverse(2,1)   &
                   +dy*dy * Sigma2_Inverse(2,2) ) &
                  / h**2
            Write (6,*)' u_j=',u_j
            K_u_j = N_norm * Exp(-u_j/Two)
            Write (6,*)' K_u_j=',K_u_j
            tmp=tmp+K_u_j
         End Do

         pxy(i,l)=tmp/DBLE(nPoints)

      End Do
   End If
End Do

Call RecPrt('pxy',' ',pxy,nPoints,nInter)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assemble the Mutual Information, MI(l), for each variable l

Do l = 1, nInter
   tmp=0.0D0
   Do i = 1, nPoints
      tmp = tmp + pxy(i,l) * log10( pxy(i,l) / (px(i,l)*py(i)) )
   End Do
   MI(l)=tmp
End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Write (*,*) 'Mutual Information'
Write (*,*) '=================='
Do i = 1, nInter
   Write (6,'(i4,E10.3)') i, MI(i)
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
Call mma_deallocate(Mean_univariate)

END SUBROUTINE pgek
