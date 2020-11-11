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
!#define _DEBUGPRINT_
SUBROUTINE pgek()
use kriging_mod, only: x, y, nPoints, nInter, Index_PGEK, nInter_Eff
Implicit None
#include "real.fh"
#include "stdalloc.fh"
Real*8, Allocatable:: Mean_univariate(:)
Real*8, Allocatable:: Variance_univariate(:), Variance_bivariate(:)
!#define _MI_SORT_
#ifdef _MI_SORT_
Integer k
#endif
Integer i, j, l
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
Real*8 d, h, u_j, K_u_j, Det_Sigma2, Sigma_2, N_norm

#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'Experimental code to compute the Mutual information between the lth coordinate'
Write (6,*) 'and the energy, denoted MI(l).'
Write (6,*)
Write (6,*) '# of sample points, nPoints=               ',nPoints
Write (6,*) '# of dimensions of the coordinates, nInter=',nInter
Call RecPrt('Energies (nPoints): y',' ',y,1,nPoints)
Call RecPrt('Coordinates (nInter x nPoints): x',' ',x,nInter,nPoints)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    Allocate memory

Call mma_allocate(Mean_univariate,nInter+1,Label='Mean_univariate')
Call mma_allocate(Variance_univariate,nInter+1,Label='Variance_univariate')
Call mma_allocate(Variance_bivariate,nInter,Label='Variance_bivariate')
Call mma_allocate(MI,nInter,Label='MI')
Call mma_allocate(px,nInter,nPoints,Label='px')
Call mma_allocate(py,nPoints,Label='py')
Call mma_allocate(pxy,nInter,nPoints,Label='pxy')

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the universal kernel density estimators px, py, and pxy
! also known as Parzen window technique using as a kernel function a normal distribution with
! an exponent of alpha=1/(2 h^2 Sigma^2), in which h is the Gaussian bandwith and Sigma^2 is the
! variance of the random variable under consideration -- Sigma is the standard deviation.

! 1) the univariate cases

! Note that for the bandwidth we are using the equation suggested by Moon, Rajagoplan and Lall
! rather than the one suggested by Liming, Qiu, Gao, Jiang, and Yang. In particular, the latter is
! identical in the uni- and bi-variate case while the former differ.
!
!h = (4/((2*d + 1)*DBLE(nPoints)))**(1/(4+d)) ! the Gaussian bandwidth   d=1,2
!h = (4/((d_h + 2)*DBLE(nPoints)))**(1/(4+d_h)) ! the Gaussian bandwidth   d_h=2
!
d=One
h = (Four/((Two*d + One)*DBLE(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
Fact = One /  ( (Two*Pi)**(d/Two) * h**d )  ! Part of the normalization factor
#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'univariate probabilities'
Write (6,*)
Write (6,*) '1/((2pi)^(d/2) h^d)=',Fact
Write (6,*) 'The Gaussian bandwidth, h=',h
#endif


!  a) Compute px(l,i)

Do l = 1, nInter
   Sigma_2=Variance_Univariate(l)
!  Write (6,*) 'l, Sigma_2=',l, Sigma_2
!
!  Take special care of coordinates which might not change due to symmetry
!  or for some other reason. Then the kernel function effectively is a Dirac delta function.
!
   If (Sigma_2<1.0D-10) Then
      px(l,:)=One/DBLE(nPoints)
   Else
      N_norm= (Fact / SQRT(Sigma_2))
      Do i = 1, nPoints
!        Write (6,*) 'i=',i

         tmp=Zero
         Do j = 1, nPoints
            u_j = (x(l,i)-x(l,j))**2 / (h**2 * Sigma_2)
            K_u_j = N_norm * EXP(-u_j/Two)
!           Write (6,*) 'u_j, K_u_j=',u_j, K_u_j
            tmp=tmp+K_u_j
         End Do

         px(l,i)=tmp/DBLE(nPoints)

      End Do
   End If
End Do

!  b) Compute py(i)

Sigma_2=Variance_Univariate(nInter+1)
N_norm= (Fact / SQRT(Sigma_2))
!Write (6,*) 'Sigma_2=',Sigma_2
Do i = 1, nPoints
!    Write (6,*) 'i=',i

   tmp=Zero
   Do j = 1, nPoints
      u_j = (y(i)-y(j))**2 / (h**2 * Sigma_2)
      K_u_j = N_norm * EXP(-u_j/Two)
!      Write (6,*) 'u_j, K_u_j=',u_j, K_u_j
      tmp=tmp+K_u_j
   End Do

   py(i)=tmp/DBLE(nPoints)

End Do
#ifdef _DEBUGPRINT_
Call RecPrt('Probability px(i,l)',' ',px,nInter,nPoints)
Call RecPrt('Probability py(l)',' ',py,nPoints,1)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2) the bivariate case, pxy(l,i)

d=Two
h = (Four/((Two*d + One)*DBLE(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
Fact = One / ( (Two*Pi)**(d/Two) * h**d )
#ifdef _DEBUGPRINT_
Write (6,*)
Write (6,*) 'bivariate probabilities'
Write (6,*)
Write (6,*) 'd=',d
Write (6,*) 'Fact=',Fact
Write (6,*) 'The Gaussian bandwidth, h=',h
#endif

Do l = 1, nInter
!  Assemble the 2 x 2 covariance matrix
   Sigma2(1,1)=Variance_univariate(l)  ! p(x_l^2)
   Sigma2(1,2)=Variance_bivariate(l)   ! p(x_ly)
   Sigma2(2,1)=Sigma2(1,2)
   Sigma2(2,2)=Variance_univariate(nInter+1) ! p(y^2)
   Det_Sigma2=Sigma2(1,1)*Sigma2(2,2)-Sigma2(1,2)*Sigma2(2,1)  ! The determinant of the 2 x 2 sigma matrix
#ifdef _DEBUGPRINT_
   Call RecPrt('Sigma2','(2E12.4)',Sigma2,2,2)
   Write (6,'(A,E12.4)') 'Det(Sigma2)=           ',Det_Sigma2
   Write (6,'(A,E12.4)') 'Sqrt(ABS(Det(Sigma2)))=',SQRT(ABS(Det_Sigma2))
#endif
   If ( Variance_univariate(l)<1.0D-10 .or. ABS(Det_Sigma2)<1.0D-20 ) Then
      pxy(l,:)=One/DBLE(nPoints)
   Else
      Sigma2_Inverse(1,1)= Sigma2(2,2)/Det_Sigma2
      Sigma2_Inverse(1,2)=-Sigma2(1,2)/Det_Sigma2
      Sigma2_Inverse(2,1)=-Sigma2(2,1)/Det_Sigma2
      Sigma2_Inverse(2,2)= Sigma2(1,1)/Det_Sigma2
!     Call RecPrt('Sigma2_Inverse','(2E10.2)',Sigma2_Inverse,2,2)
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
!           Write (6,*)' u_j=',u_j
            K_u_j = N_norm * Exp(-u_j/Two)
!           Write (6,*)' K_u_j=',K_u_j
            tmp=tmp+K_u_j
         End Do

         pxy(l,i)=tmp/DBLE(nPoints)

      End Do
   End If
!  Write (*,*) 'pxy(l,:)=',pxy(l,:)
End Do

#ifdef _DEBUGPRINT_
Call RecPrt('Probability pxy(l,i)',' ',pxy,nInter,nPoints)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assemble the Mutual Information, MI(l), for each variable l, according to equation 23 of
! Chen et al. Note that this is correct only if the sample points to a large extent are
! equidistant -- which they are not.

Do l = 1, nInter
   tmp=Zero
   If (Variance_univariate(l)>1.0D-10) Then
      Do i = 1, nPoints
         tmp = tmp + pxy(l,i) * log10( pxy(l,i) / (  px(l,i) * py(i)) )
      End Do
   End If
   MI(l)=tmp
End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _DEBUGPRINT_
Write (6,*) 'Mutual Information'
Write (6,*) '=================='
Do i = 1, nInter
   Write (6,'(i4,E10.3)') i, MI(i)
End Do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nInter_Eff=0
Do i = 1, nInter
#ifdef _MI_SORT_
   k = 0
   tmp=-One
   Do j = 1, nInter
      If (MI(j)>tmp) Then
         k=j
         tmp=MI(j)
      End If
   End Do
   If (k==0) Then
      Write (6,*) 'PGEK: k==0'
      Call Abend()
   Else
      If (MI(k)>1.0D-10) nInter_Eff=nInter_Eff+1
      MI(k)=-Two
      Index_PGEK(i)=k
   End If
#else
!  If (Abs(Variance_bivariate(i))>5.0D-15) Then
   If (Abs(Variance_bivariate(i))>1.0D-14) Then
      nInter_eff=nInter_eff+1
      Index_PGEK(nInter_eff)=i
   End If
#endif
End Do
!
! code debugging with full lists and forward and backward order.
!
!nInter_Eff = nInter
!Do i = 1, nInter
!   Index_PGEK(i)=nInter-i+1
!End Do
!nInter_Eff = nInter
!Do i = 1, nInter
!   Index_PGEK(i)=i
!End Do
#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
Write (6,*) 'nInter_eff=',nInter_eff
Write (6,*) 'Index_PGEK=',(Index_PGEK(i),i=1,nInter_eff)
#endif
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

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real*8 Function p_x(xl,l)
Real*8 :: xl
Integer l

d=One
h = (Four/((Two*d + One)*DBLE(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
Fact = One /  ( (Two*Pi)**(d/Two) * h**d )  ! Part of the normalization factor

Sigma_2=Variance_Univariate(l)

! We need a case of explicit zero variance (within machine accuracy) for the cases for symmetry breaking coordinates.

N_norm= (Fact / SQRT(Sigma_2))
p_x=Zero
Do j = 1, nPoints
   u_j = (xl-x(l,j))**2 / (h**2 * Sigma_2)
   K_u_j = N_norm * EXP(-u_j/Two)
   p_x=p_x+K_u_j
End Do
p_x = p_x / DBLE(nPoints)
End Function p_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real*8 Function p_y(yl)
Real*8 :: yl

d=One
h = (Four/((Two*d + One)*DBLE(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
Fact = One /  ( (Two*Pi)**(d/Two) * h**d )  ! Part of the normalization factor

Sigma_2=Variance_Univariate(nPoints+1)
N_norm= (Fact / SQRT(Sigma_2))
p_y=Zero
Do j = 1, nPoints
   u_j = (yl-y(j))**2 / (h**2 * Sigma_2)
   K_u_j = N_norm * EXP(-u_j/Two)
   p_y=p_y+K_u_j
End Do
p_y = p_y / DBLE(nPoints)

End Function p_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real*8 Function p_xy(xl,yl,l)
Real*8 :: xl, yl
Integer l

d=Two
h = (Four/((Two*d + One)*DBLE(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
Fact = One / ( (Two*Pi)**(d/Two) * h**d )

! Assemble the 2 x 2 covariance matrix
Sigma2(1,1)=Variance_univariate(l)  ! p(x_l^2)
Sigma2(1,2)=Variance_bivariate(l)   ! p(x_ly)
Sigma2(2,1)=Sigma2(1,2)
Sigma2(2,2)=Variance_univariate(nInter+1) ! p(y^2)
Det_Sigma2=Sigma2(1,1)*Sigma2(2,2)-Sigma2(1,2)*Sigma2(2,1)  ! The determinant of the 2 x 2 sigma matrix
Sigma2_Inverse(1,1)= Sigma2(2,2)/Det_Sigma2
Sigma2_Inverse(1,2)=-Sigma2(1,2)/Det_Sigma2
Sigma2_Inverse(2,1)=-Sigma2(2,1)/Det_Sigma2
Sigma2_Inverse(2,2)= Sigma2(1,1)/Det_Sigma2

N_norm = (Fact / SQRT(Det_Sigma2))
p_xy=Zero
Do j = 1, nPoints
   dx = xl-x(l,j)
   dy = yl-y(j)
   u_j = ( dx*dx * Sigma2_Inverse(1,1)   &
           +dx*dy * Sigma2_Inverse(1,2)   &
           +dy*dx * Sigma2_Inverse(2,1)   &
           +dy*dy * Sigma2_Inverse(2,2) ) &
       / h**2
   K_u_j = N_norm * Exp(-u_j/Two)
   p_xy=p_xy+K_u_j
End Do
p_xy = p_xy / DBLE(nPoints)

End Function p_xy

END SUBROUTINE pgek
