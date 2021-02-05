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

subroutine pgek()

use kriging_mod,only:x,y,nPoints,nInter,Index_PGEK,nInter_Eff

implicit none
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: Mean_univariate(:)
real*8, allocatable :: Variance_univariate(:), Variance_bivariate(:)
!#define _MI_SORT_
#ifdef _MI_SORT_
integer k
#endif
integer i, j, l
real*8 tmp, dx, dy, Fact
! Mutual information array
!real*8 MI(nInter)
real*8, allocatable :: MI(:)
! universal kernal density estimators.
!real*8 px(nPoints,nInter)
!real*8 py(nPoints)
!real*8 pxy(nPoints,nInter)
real*8, allocatable :: px(:,:), py(:), pxy(:,:)
! the 2 x 2 covariance matrix (used to compute the transformation variable in the bivariate case)
real*8 Sigma2(2,2), Sigma2_Inverse(2,2)
real*8 d, h, u_j, K_u_j, Det_Sigma2, Sigma_2, N_norm

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'Experimental code to compute the Mutual information between the lth coordinate'
write(6,*) 'and the energy, denoted MI(l).'
write(6,*)
write(6,*) '# of sample points, nPoints=               ',nPoints
write(6,*) '# of dimensions of the coordinates, nInter=',nInter
call RecPrt('Energies (nPoints): y',' ',y,1,nPoints)
call RecPrt('Coordinates (nInter x nPoints): x',' ',x,nInter,nPoints)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allocate memory

call mma_allocate(Mean_univariate,nInter+1,Label='Mean_univariate')
call mma_allocate(Variance_univariate,nInter+1,Label='Variance_univariate')
call mma_allocate(Variance_bivariate,nInter,Label='Variance_bivariate')
call mma_allocate(MI,nInter,Label='MI')
call mma_allocate(px,nInter,nPoints,Label='px')
call mma_allocate(py,nPoints,Label='py')
call mma_allocate(pxy,nInter,nPoints,Label='pxy')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute the mean for all variables - coordinate components and the energy
! univariate --  x and y
! Compute the variance of the univariate and the bivariate

! univariate means
do i=1,nInter
  Mean_Univariate(i) = Zero
  do j=1,nPoints
    Mean_Univariate(i) = Mean_Univariate(i)+x(i,j)
  end do
end do
Mean_Univariate(nInter+1) = Zero
do j=1,nPoints
  Mean_Univariate(nInter+1) = Mean_Univariate(nInter+1)+y(j)
end do
Mean_Univariate(:) = Mean_Univariate(:)/dble(nPoints)

! uni- and bivariate variances
do i=1,nInter
  Variance_Univariate(i) = Zero
  Variance_bivariate(i) = Zero
  do j=1,nPoints
    Variance_Univariate(i) = Variance_Univariate(i)+(x(i,j)-Mean_univariate(i))**2
    Variance_bivariate(i) = Variance_bivariate(i)+(x(i,j)-Mean_univariate(i))*(y(j)-Mean_univariate(nInter+1))
  end do
end do
Variance_Univariate(nInter+1) = Zero
do j=1,nPoints
  Variance_Univariate(nInter+1) = Variance_Univariate(nInter+1)+(y(j)-Mean_univariate(nInter+1))**2
end do
Variance_Univariate(:) = Variance_Univariate(:)/dble(nPoints)
Variance_Bivariate(:) = Variance_Bivariate(:)/dble(nPoints)

#ifdef _DEBUGPRINT_
call RecPrt('Mean - x,y ',' ',Mean_Univariate(:),nInter+1,1)
call RecPrt('Variance - x^2,y^2 ',' ',Variance_Univariate(:),nInter+1,1)
call RecPrt('Variance - xy ',' ',Variance_Bivariate(:),nInter,1)
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

d = One
h = (Four/((Two*d+One)*dble(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
Fact = One/((Two*Pi)**(d/Two)*h**d)  ! Part of the normalization factor
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'univariate probabilities'
write(6,*)
write(6,*) '1/((2pi)^(d/2) h^d)=',Fact
write(6,*) 'The Gaussian bandwidth, h=',h
#endif

! a) Compute px(l,i)

do l=1,nInter
  Sigma_2 = Variance_Univariate(l)
  !write (6,*) 'l, Sigma_2=',l, Sigma_2
  !
  ! Take special care of coordinates which might not change due to symmetry
  ! or for some other reason. Then the kernel function effectively is a Dirac delta function.

  if (Sigma_2 < 1.0d-10) then
    px(l,:) = One/dble(nPoints)
  else
    N_norm = (Fact/sqrt(Sigma_2))
    do i=1,nPoints
      !write (6,*) 'i=',i

      tmp = Zero
      do j=1,nPoints
        u_j = (x(l,i)-x(l,j))**2/(h**2*Sigma_2)
        K_u_j = N_norm*exp(-u_j/Two)
        !write (6,*) 'u_j, K_u_j=',u_j, K_u_j
        tmp = tmp+K_u_j
      end do

      px(l,i) = tmp/dble(nPoints)

    end do
  end if
end do

! b) Compute py(i)

Sigma_2 = Variance_Univariate(nInter+1)
N_norm = (Fact/sqrt(Sigma_2))
!write (6,*) 'Sigma_2=',Sigma_2
do i=1,nPoints
  !write (6,*) 'i=',i

  tmp = Zero
  do j=1,nPoints
    u_j = (y(i)-y(j))**2/(h**2*Sigma_2)
    K_u_j = N_norm*exp(-u_j/Two)
    !write (6,*) 'u_j, K_u_j=',u_j, K_u_j
    tmp = tmp+K_u_j
  end do

  py(i) = tmp/dble(nPoints)

end do
#ifdef _DEBUGPRINT_
call RecPrt('Probability px(i,l)',' ',px,nInter,nPoints)
call RecPrt('Probability py(l)',' ',py,nPoints,1)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2) the bivariate case, pxy(l,i)

d = Two
h = (Four/((Two*d+One)*dble(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
Fact = One/((Two*Pi)**(d/Two)*h**d)
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'bivariate probabilities'
write(6,*)
write(6,*) 'd=',d
write(6,*) 'Fact=',Fact
write(6,*) 'The Gaussian bandwidth, h=',h
#endif

do l=1,nInter
  ! Assemble the 2 x 2 covariance matrix
  Sigma2(1,1) = Variance_univariate(l)  ! p(x_l^2)
  Sigma2(1,2) = Variance_bivariate(l)   ! p(x_ly)
  Sigma2(2,1) = Sigma2(1,2)
  Sigma2(2,2) = Variance_univariate(nInter+1) ! p(y^2)
  Det_Sigma2 = Sigma2(1,1)*Sigma2(2,2)-Sigma2(1,2)*Sigma2(2,1)  ! The determinant of the 2 x 2 sigma matrix
# ifdef _DEBUGPRINT_
  call RecPrt('Sigma2','(2E12.4)',Sigma2,2,2)
  write(6,'(A,E12.4)') 'Det(Sigma2)=           ',Det_Sigma2
  write(6,'(A,E12.4)') 'Sqrt(ABS(Det(Sigma2)))=',sqrt(abs(Det_Sigma2))
# endif
  if (Variance_univariate(l) < 1.0d-10 .or. abs(Det_Sigma2) < 1.0d-20) then
    pxy(l,:) = One/dble(nPoints)
  else
    Sigma2_Inverse(1,1) = Sigma2(2,2)/Det_Sigma2
    Sigma2_Inverse(1,2) = -Sigma2(1,2)/Det_Sigma2
    Sigma2_Inverse(2,1) = -Sigma2(2,1)/Det_Sigma2
    Sigma2_Inverse(2,2) = Sigma2(1,1)/Det_Sigma2
    !call RecPrt('Sigma2_Inverse','(2E10.2)',Sigma2_Inverse,2,2)
    N_norm = (Fact/sqrt(Det_Sigma2))
    do i=1,nPoints

      tmp = Zero
      do j=1,nPoints
        dx = x(l,i)-x(l,j)
        dy = y(i)-y(j)
        u_j = (dx*dx*Sigma2_Inverse(1,1)+dx*dy*Sigma2_Inverse(1,2)+dy*dx*Sigma2_Inverse(2,1)+dy*dy*Sigma2_Inverse(2,2))/h**2
        !write (6,*)' u_j=',u_j
        K_u_j = N_norm*exp(-u_j/Two)
        !write (6,*)' K_u_j=',K_u_j
        tmp = tmp+K_u_j
      end do

      pxy(l,i) = tmp/dble(nPoints)

    end do
  end if
  !write (6,*) 'pxy(l,:)=',pxy(l,:)
end do

#ifdef _DEBUGPRINT_
call RecPrt('Probability pxy(l,i)',' ',pxy,nInter,nPoints)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assemble the Mutual Information, MI(l), for each variable l, according to equation 23 of
! Chen et al. Note that this is correct only if the sample points to a large extent are
! equidistant -- which they are not.

do l=1,nInter
  tmp = Zero
  if (Variance_univariate(l) > 1.0d-10) then
    do i=1,nPoints
      tmp = tmp+pxy(l,i)*log10(pxy(l,i)/(px(l,i)*py(i)))
    end do
  end if
  MI(l) = tmp
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _DEBUGPRINT_
write(6,*) 'Mutual Information'
write(6,*) '=================='
do i=1,nInter
  write(6,'(i4,E10.3)') i,MI(i)
end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nInter_Eff = 0
do i=1,nInter
# ifdef _MI_SORT_
  k = 0
  tmp = -One
  do j=1,nInter
    if (MI(j) > tmp) then
      k = j
      tmp = MI(j)
    end if
  end do
  if (k == 0) then
    write(6,*) 'PGEK: k==0'
    call Abend()
  else
    if (MI(k) > 1.0d-10) nInter_Eff = nInter_Eff+1
    MI(k) = -Two
    Index_PGEK(i) = k
  end if
# else
  ! If (Abs(Variance_bivariate(i))>5.0D-15) Then
  if (abs(Variance_bivariate(i)) > 1.0d-14) then
    nInter_eff = nInter_eff+1
    Index_PGEK(nInter_eff) = i
  end if
# endif
end do

! code debugging with full lists and forward and backward order.

!nInter_Eff = nInter
!do i=1, nInter
! Index_PGEK(i) = nInter-i+1
!end do
!nInter_Eff = nInter
!do i=1, nInter
! Index_PGEK(i) = i
!end do
#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(6,*) 'nInter_eff=',nInter_eff
write(6,*) 'Index_PGEK=',(Index_PGEK(i),i=1,nInter_eff)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Deallocate memory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_deallocate(MI)
call mma_deallocate(px)
call mma_deallocate(py)
call mma_deallocate(pxy)
call mma_deallocate(Variance_univariate)
call mma_deallocate(Variance_bivariate)
call mma_deallocate(Mean_univariate)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function p_x(xl,l)

  real*8 :: xl
  integer :: j,l

  d = One
  h = (Four/((Two*d+One)*dble(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
  Fact = One/((Two*Pi)**(d/Two)*h**d)  ! Part of the normalization factor

  Sigma_2 = Variance_Univariate(l)

  ! We need a case of explicit zero variance (within machine accuracy) for the cases for symmetry breaking coordinates.

  N_norm = (Fact/sqrt(Sigma_2))
  p_x = Zero
  do j=1,nPoints
    u_j = (xl-x(l,j))**2/(h**2*Sigma_2)
    K_u_j = N_norm*exp(-u_j/Two)
    p_x = p_x+K_u_j
  end do
  p_x = p_x/dble(nPoints)
end function p_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function p_y(yl)

  real*8 :: yl
  integer :: j

  d = One
  h = (Four/((Two*d+One)*dble(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
  Fact = One/((Two*Pi)**(d/Two)*h**d)  ! Part of the normalization factor

  Sigma_2 = Variance_Univariate(nPoints+1)
  N_norm = (Fact/sqrt(Sigma_2))
  p_y = Zero
  do j=1,nPoints
    u_j = (yl-y(j))**2/(h**2*Sigma_2)
    K_u_j = N_norm*exp(-u_j/Two)
    p_y = p_y+K_u_j
  end do
  p_y = p_y/dble(nPoints)

end function p_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function p_xy(xl,yl,l)

  real*8 :: xl,yl
  integer :: j,l

  d = Two
  h = (Four/((Two*d+One)*dble(nPoints)))**(One/(Four+d)) ! the Gaussian bandwidth
  Fact = One/((Two*Pi)**(d/Two)*h**d)

  ! Assemble the 2 x 2 covariance matrix
  Sigma2(1,1) = Variance_univariate(l)  ! p(x_l^2)
  Sigma2(1,2) = Variance_bivariate(l)   ! p(x_ly)
  Sigma2(2,1) = Sigma2(1,2)
  Sigma2(2,2) = Variance_univariate(nInter+1) ! p(y^2)
  Det_Sigma2 = Sigma2(1,1)*Sigma2(2,2)-Sigma2(1,2)*Sigma2(2,1)  ! The determinant of the 2 x 2 sigma matrix
  Sigma2_Inverse(1,1) = Sigma2(2,2)/Det_Sigma2
  Sigma2_Inverse(1,2) = -Sigma2(1,2)/Det_Sigma2
  Sigma2_Inverse(2,1) = -Sigma2(2,1)/Det_Sigma2
  Sigma2_Inverse(2,2) = Sigma2(1,1)/Det_Sigma2

  N_norm = (Fact/sqrt(Det_Sigma2))
  p_xy = Zero
  do j=1,nPoints
    dx = xl-x(l,j)
    dy = yl-y(j)
    u_j = (dx*dx*Sigma2_Inverse(1,1)+dx*dy*Sigma2_Inverse(1,2)+dy*dx*Sigma2_Inverse(2,1)+dy*dy*Sigma2_Inverse(2,2))/h**2
    K_u_j = N_norm*exp(-u_j/Two)
    p_xy = p_xy+K_u_j
  end do
  p_xy = p_xy/dble(nPoints)

end function p_xy

end subroutine pgek
