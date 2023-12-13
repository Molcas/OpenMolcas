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
! Copyright (C) 2013,2014, Ignacio Fdez. Galvan                        *
!***********************************************************************
!> @file
!> @brief
!>   Subroutines for superposition of molecules.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Subroutines for superposition of molecules.
!> \cite The2005-ACSA-61-478
!> \cite Liu2010-JCC-31-1561
!***********************************************************************

!***********************************************************************
!  Get_RMSD_w
!
!> @brief
!>   Compute the optimal weighted RMSD between two structures.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Compute the optimal RMSD between two structures, with weights.
!> Coordinates are not changed, but the RMSD corresponds to the aligned structures.
!>
!> @param[in,out] x    Cartesian coordinates of the first structure
!> @param[in]     y    Cartesian coordinates of the second structure
!> @param[in]     w    Weights for each atom
!> @param[in]     nAt  Number of atoms in the structures
!> @param[out]    RMSD Minimum weighted RMSD between the two structures
!>
!> @see ::Get_RMSD
!***********************************************************************
subroutine Get_RMSD_w(x,y,w,nAt,RMSD)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: x(3,nAt)
real(kind=wp), intent(in) :: y(3,nAt), w(nAt)
real(kind=wp), intent(out) :: RMSD

call get_rotation(x,y,w,nAt,RMSD,.false.)

end subroutine Get_RMSD_w

!***********************************************************************
!  Get_RMSD
!
!> @brief
!>   Compute the optimal RMSD between two structures.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Compute the optimal RMSD between two structures, with equal weights.
!> Coordinates are not changed, but the RMSD corresponds to the aligned structures.
!>
!> @param[in,out] x    Cartesian coordinates of the first structure
!> @param[in]     y    Cartesian coordinates of the second structure
!> @param[in]     nAt  Number of atoms in the structures
!> @param[out]    RMSD Minimum RMSD between the two structures
!>
!> @see ::Get_RMSD_w
!***********************************************************************
subroutine Get_RMSD(x,y,nAt,RMSD)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: x(3,nAt)
real(kind=wp), intent(in) :: y(3,nAt)
real(kind=wp), intent(out) :: RMSD
real(kind=wp), allocatable :: w(:)

call mma_allocate(w,nAt)
w(:) = One
call Get_RMSD_w(x,y,w,nAt,RMSD)
call mma_deallocate(w)

end subroutine Get_RMSD

!***********************************************************************
!  Superpose_w
!
!> @brief
!>   Superpose two structures.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Superpose two structures, with weights.
!> The coordinates of the first structure are modified.
!>
!> @param[in,out]  x    Cartesian coordinates of the first structure
!> @param[in]      y    Cartesian coordinates of the second structure
!> @param[in]      w    Weights for each atom
!> @param[in]      nAt  Number of atoms in the structures
!> @param[out]     RMSD Weighted RMSD between the two final structures
!> @param[out]     RMax Maximum atom-wise (weighted) distance between the two structures
!>
!> @see ::Superpose
!***********************************************************************
subroutine Superpose_w(x,y,w,nAt,RMSD,RMax)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: x(3,nAt)
real(kind=wp), intent(in) :: y(3,nAt), w(nAt)
real(kind=wp), intent(out) :: RMSD, RMax
integer(kind=iwp) :: iAt
real(kind=wp) :: r2

call get_rotation(x,y,w,nAt,RMSD,.true.)
! Calculate the maximum weighted distance between atoms
! (not sure what meaning it has with w!=1)
RMax = Zero
do iAt=1,nAt
  r2 = (x(1,iAt)-y(1,iAt))**2+(x(2,iAt)-y(2,iAt))**2+(x(3,iAt)-y(3,iAt))**2
  RMax = max(RMax,w(iAt)*r2)
end do
RMax = sqrt(RMax)

end subroutine Superpose_w

!***********************************************************************
!  Superpose_w
!
!> @brief
!>   Superpose two structures.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Superpose two structures, with equal weights.
!> The coordinates of the first structure are modified.
!>
!> @param[in,out]  x    Cartesian coordinates of the first structure
!> @param[in]      y    Cartesian coordinates of the second structure
!> @param[in]      nAt  Number of atoms in the structures
!> @param[out]     RMSD RMSD between the two final structures
!> @param[out]     RMax Maximum atom-wise distance between the two structures
!>
!> @see ::Superpose_w
!***********************************************************************
subroutine Superpose(x,y,nAt,RMSD,RMax)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: x(3,nAt)
real(kind=wp), intent(in) :: y(3,nAt)
real(kind=wp), intent(out) :: RMSD, RMax
real(kind=wp), allocatable :: w(:)

call mma_allocate(w,nAt)
w(:) = One
call Superpose_w(x,y,w,nAt,RMSD,RMax)
call mma_deallocate(w)

end subroutine Superpose

!***********************************************************************
!> @name Internal
!>
!> These procedures are probably not very useful outside the file.
!> @{
!***********************************************************************
!  get_rotation
!
!> @brief
!>   Compute the RMSD and superpose two structures.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Calculate the minimum weighted RMSD between two structures and, optionally,
!> rotate and translate the first one to be superposed onto the other.
!>
!> @param[in,out] x      Cartesian coordinates of the first structure
!> @param[in]     y      Cartesian coordinates of the second structure
!> @param[in]     w      Weights for each atom
!> @param[in]     nAt    Number of atoms in the structures
!> @param[out]    RMSD   Minimum weighted RMSD between the two final structures
!> @param[in]     rotate Flag for changing the coordinates or not
!***********************************************************************
subroutine get_rotation(x,y,w,nAt,RMSD,rotate)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
! if USE_QCD is defined, use the real QCD method, otherwise use conventional
! method to locate the largest eigenvalue (may be more robust in some cases)
!#define USE_QCD
#ifdef USE_QCD
use Constants, only: Half
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: x(3,nAt)
real(kind=wp), intent(in) :: y(3,nAt), w(nAt)
real(kind=wp), intent(out) :: RMSD
logical(kind=iwp), intent(in) :: rotate
integer(kind=iwp) :: i, iAt
real(kind=wp) :: c_x(3), c_y(3), Gx, Gy, Kxy(4,4), lambda, Mxy(3,3), q(4), wTot
#ifdef USE_QCD
real(kind=wp) :: Poly(5)
#else
real(kind=wp) :: wTmp(100)
#endif
real(kind=wp), allocatable :: xCen(:,:), yCen(:,:)
real(kind=wp), external :: DDot_, inner_prod

call mma_allocate(xCen,3,nAt)
call mma_allocate(yCen,3,nAt)
! Center the two structures and calculate their inner products
call center_mol(x,w,nAt,c_x,xCen)
call center_mol(y,w,nAt,c_y,yCen)
Gx = inner_prod(xCen,w,nAt)
Gy = inner_prod(yCen,w,nAt)
call inner_mat(xCen,yCen,w,nAt,Mxy)
! Find the largest eigenvalue (always lower or equal to (Gx+Gy)/2)
#ifdef USE_QCD
call build_polynomial(Mxy,Poly)
lambda = Half*(Gx+Gy)
call find_lambda(Poly,lambda)
#else
call build_K_matrix(Mxy,Kxy)
i = 0
call dsyev_('N','U',4,Kxy,4,q,wTmp,100,i)
lambda = q(4)
#endif
! Calculate the optimized RMSD value
wTot = sum(w)
RMSD = sqrt(abs(Gx+Gy-Two*lambda)/wTot)
! The rotation is only performed if an actual superposition is requested
if (rotate) then
  ! The rotation quaternion is obtained as the eigenvector of K corresponding
  ! to the eigenvalue lambda.
  call build_K_matrix(Mxy,Kxy)
  call get_eigenvector(Kxy,lambda,q)
  ! Apply the rotation to the first structure. Before rotation, the quaternion
  ! must be normalized and the rotation inverted (change of sign in the 1st component)
  q(:) = q/sqrt(DDot_(4,q,1,q,1))
  q(1) = -q(1)
  call apply_rotation(xCen,nAt,q)
  ! Translate the structure to match the second (reference) one
  do iAt=1,nAt
    xCen(:,iAt) = xCen(:,iAt)+c_y
  end do
  ! Copy the aligned structures back in the original matrices
  x(:,:) = xCen
end if
call mma_deallocate(xCen)
call mma_deallocate(yCen)

end subroutine get_rotation

!***********************************************************************
!  center_mol
!
!> @brief
!>   Center a structure on its weighted geometric center.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Center a structure on its weighted geometric center.
!>
!> @param[in]  x    Cartesian coordinates of the structure
!> @param[in]  w    Weights for each atom
!> @param[in]  nAt  Number of atoms in the structures
!> @param[out] c    Weighted center
!> @param[out] xCen Centered Cartesian coordinates
!***********************************************************************
subroutine center_mol(x,w,nAt,c,xCen)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(in) :: x(3,nAt), w(nAt)
real(kind=wp), intent(out) :: c(3), xCen(3,nAt)
integer(kind=iwp) :: i, iAt
real(kind=wp) :: wTot

! Compute the total weight
wTot = sum(w)
do i=1,3
  ! Calculate the center in each dimension (x, y, z)
  c(i) = Zero
  do iAt=1,nAt
    c(i) = c(i)+w(iAt)*x(i,iAt)
  end do
  c(i) = c(i)/wTot
  ! Center the structure in each dimension
  xCen(i,:) = x(i,:)-c(i)
end do

end subroutine center_mol

!***********************************************************************
!  inner_prod
!
!> @brief
!>   Calculate the inner product of a single structure.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Returns the inner product of the weighted coordinates of a single structure:
!>  \f[P = \sum_{i=1,N} w_i ( x_i^2 + y_i^2 + z_i^2 ) \f]
!>
!> @param[in] x   Cartesian coordinates of the structure
!> @param[in] w   Weights for each atom
!> @param[in] nAt Number of atoms in the structures
!>
!> @return The inner product \f$ P \f$
!***********************************************************************
function inner_prod(x,w,nAt)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: inner_prod
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(in) :: x(3,nAt), w(nAt)
integer(kind=iwp) :: iAt
real(kind=wp) :: r

r = Zero
do iAt=1,nAt
  r = r+w(iAt)*(x(1,iAt)**2+x(2,iAt)**2+x(3,iAt)**2)
end do
inner_prod = r

end function inner_prod

!***********************************************************************
!  inner_mat
!
!> @brief
!>   Calculate the inner product matrix of two structures.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Calculate the inner product matrix of the weighted coordinates of two structures,
!> \f$ A, B \f$:
!>   \f[ M_{ij} = \sum_{k=1,N} w_i (A_i)_k (B_j)_k \quad i,j \in \{x, y, z\} \f]
!>
!> @param[in]  x   Cartesian coordinates of the first structure
!> @param[in]  y   Cartesian coordinates of the second structure
!> @param[in]  w   Weights for each atom
!> @param[in]  nAt Number of atoms in the structures
!> @param[out] M   Inner product matrix
!***********************************************************************
subroutine inner_mat(x,y,w,nAt,M)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(in) :: x(3,nAt), y(3,nAt), w(nAt)
real(kind=wp), intent(out) :: M(3,3)
integer(kind=iwp) :: i, iAt

M(:,:) = Zero
do iAt=1,nAt
  do i=1,3
    M(:,i) = M(:,i)+w(iAt)*x(:,iAt)*y(i,iAt)
  end do
end do

end subroutine inner_mat

!***********************************************************************
!  build_K_matrix
!
!> @brief
!>   Build the key matrix \f$ K \f$ from the \f$ M \f$ matrix.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Build the key matrix \f$ K \f$ from the inner product matrix \f$ M \f$ matrix.
!>
!> @param[in]  M Inner product matrix
!> @param[out] K Key matrix \f$ K \f$
!***********************************************************************
subroutine build_K_matrix(M,K)

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: M(3,3)
real(kind=wp), intent(out) :: K(4,4)
integer(kind=iwp) :: i, j

! Compute the unique elements
K(1,1) = M(1,1)+M(2,2)+M(3,3)
K(1,2) = M(2,3)-M(3,2)
K(1,3) = M(3,1)-M(1,3)
K(1,4) = M(1,2)-M(2,1)
K(2,2) = M(1,1)-M(2,2)-M(3,3)
K(2,3) = M(1,2)+M(2,1)
K(2,4) = M(1,3)+M(3,1)
K(3,3) = M(2,2)-M(1,1)-M(3,3)
K(3,4) = M(2,3)+M(3,2)
K(4,4) = M(3,3)-M(1,1)-M(2,2)
! Make the matrix symmetric
do i=2,4
  do j=1,i-1
    K(i,j) = K(j,i)
  end do
end do

end subroutine build_K_matrix

!***********************************************************************
! build_polynomial
!
!> @brief
!>   Calculate the coefficients for the characteristic polynomial.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Calculate the coefficients for the characteristic polynomial of the \f$ K \f$
!> matrix directly from the elements of the \f$ M \f$ matrix.
!>
!> @param[in]  M Inner product matrix
!> @param[out] P The polynomial coefficients
!***********************************************************************
subroutine build_polynomial(M,P)

use Constants, only: Zero, One, Two, Eight
use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: M(3,3)
real(kind=wp), intent(out) :: P(5)
real(kind=wp) :: a1, a2, a3, a4, a5, a6
real(kind=wp), external :: DDot_, determinant3

! The x^4 and x^3 coefficients are constant
P(5) = One
P(4) = Zero
! The x^2 coefficient is the sum of all squared elements, times -2
P(3) = -Two*DDot_(9,M,1,M,1)
! The x^1 coefficient is the determinant of M, times -8
P(2) = -Eight*determinant3(M)
! The x^0 coefficient has a rather cumbersome expansion
a1 = (M(1,2)**2+M(1,3)**2-M(2,1)**2-M(3,1)**2)**2
a2 = (-M(1,1)**2+M(2,2)**2+M(3,3)**2+M(2,3)**2+M(3,2)**2-Two*(M(2,2)*M(3,3)-M(2,3)*M(3,2)))* &
     (-M(1,1)**2+M(2,2)**2+M(3,3)**2+M(2,3)**2+M(3,2)**2+Two*(M(2,2)*M(3,3)-M(2,3)*M(3,2)))
a3 = (-(M(1,3)+M(3,1))*(M(2,3)-M(3,2))+(M(1,2)-M(2,1))*(M(1,1)-M(2,2)-M(3,3)))* &
     (-(M(1,3)-M(3,1))*(M(2,3)+M(3,2))+(M(1,2)-M(2,1))*(M(1,1)-M(2,2)+M(3,3)))
a4 = (-(M(1,3)+M(3,1))*(M(2,3)+M(3,2))-(M(1,2)+M(2,1))*(M(1,1)+M(2,2)-M(3,3)))* &
     (-(M(1,3)-M(3,1))*(M(2,3)-M(3,2))-(M(1,2)+M(2,1))*(M(1,1)+M(2,2)+M(3,3)))
a5 = ((M(1,2)+M(2,1))*(M(2,3)+M(3,2))+(M(1,3)+M(3,1))*(M(1,1)-M(2,2)+M(3,3)))* &
     (-(M(1,2)-M(2,1))*(M(2,3)-M(3,2))+(M(1,3)+M(3,1))*(M(1,1)+M(2,2)+M(3,3)))
a6 = ((M(1,2)+M(2,1))*(M(2,3)-M(3,2))+(M(1,3)-M(3,1))*(M(1,1)-M(2,2)-M(3,3)))* &
     (-(M(1,2)-M(2,1))*(M(2,3)+M(3,2))+(M(1,3)-M(3,1))*(M(1,1)+M(2,2)-M(3,3)))
P(1) = a1+a2+a3+a4+a5+a6

end subroutine build_polynomial

!***********************************************************************
!  find_lambda
!
!> @brief
!>   Find one root of a 4th degree polynomial.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Find one real root of a 4th degree polynomial using the Newton--Raphson method.
!> On input, the variable \p r contains the initial guess for the root.
!>
!> @param[in]     P The polynomial coefficients
!> @param[in,out] r The found root (initial guess on input)
!***********************************************************************
subroutine find_lambda(P,r)

use Constants, only: Zero, Two, Ten
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: P(5)
real(kind=wp), intent(inout) :: r
integer(kind=iwp) :: j, loop
real(kind=wp) :: df, f, r_old, thrr
integer(kind=iwp), parameter :: mxloop = 100
real(kind=wp), parameter :: thr = 1.0e-11_wp

r_old = max(Two*r,Ten)
! Find the root with the Newton-Raphson method, starting from the initial value
loop = 0
thrr = thr*r
do while ((abs(r-r_old) > thrr) .and. (loop < mxloop))
  r_old = r
  ! Compute the values of the polynomial (f) and its derivative (df)
  ! at the current value of r
  f = P(5)
  df = Zero
  do j=4,1,-1
    df = df*r+f
    f = f*r+P(j)
  end do
  ! Update the value of r from the values of f and df
  ! (special case if the derivative is zero)
  if (abs(df) < thrr) then
    if (abs(f) < thr) then
      r_old = r
    else
      ! If the df is zero and r is not a root, we are in trouble,
      ! just move the point a tad
      r = r-sign(Two*thrr,f)
    end if
  else
    r = r-f/df
  end if
  thrr = thr*r
  ! Increase count to avoid neverending loop
  loop = loop+1
end do

end subroutine find_lambda

!***********************************************************************
!  get_eigenvector
!
!> @brief
!>   Find an eigenvector of the \f$ K \f$ matrix, given its eigenvalue.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Compute the eigenvector of the \f$ K \f$ matrix corresponding to one eigenvalue.
!> @side_effects
!> The matrix \p K is destroyed on output.
!>
!> @param[in,out] K The \f$ K \f$ matrix
!> @param[in]     r The eigenvalue
!> @param[out]    V The eigenvector
!***********************************************************************
subroutine get_eigenvector(K,r,V)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: K(4,4)
real(kind=wp), intent(in) :: r
real(kind=wp), intent(out) :: V(4)
integer(kind=iwp) :: i, j
real(kind=wp) :: n
real(kind=wp), parameter :: thr = 1.0e-12_wp
real(kind=wp), external :: cofactor, DDot_

! Get the matrix K-rI
K(1,1) = K(1,1)-r
K(2,2) = K(2,2)-r
K(3,3) = K(3,3)-r
K(4,4) = K(4,4)-r
! Find the eigenvector as any non-zero column of adj(K-rI)
n = Zero
do j=1,4
  if (n < thr) then
    ! Since adj(X) is the transpose of the cofactors matrix, we take a row from K
    do i=1,4
      V(i) = cofactor(K,j,i)
    end do
    ! Get the norm of the vector, if not too small, no further vector will be computed
    n = DDot_(4,V,1,V,1)
  end if
end do
! If the norm is still too small, return the identity quaternion (no rotation)
if (n < thr) V(:) = [One,Zero,Zero,Zero]

end subroutine get_eigenvector

!***********************************************************************
!  cofactor
!
!> @brief
!>   Return a cofactor from a 4&times;4 matrix.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Returns a cofactor (determinant of a submatrix) from a 4&times;4 matrix.
!>
!> @param[in] M The matrix
!> @param[in] i The row of the cofactor
!> @param[in] j The column of the cofactor
!>
!> @return The \f$ C_{ij} \f$ cofactor of the matrix \p M
!***********************************************************************
function cofactor(M,i,j)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: cofactor
real(kind=wp), intent(in) :: M(4,4)
integer(kind=iwp), intent(in) :: i, j
integer(kind=iwp) :: ii, jj
real(kind=wp) :: A(3,3), f
real(kind=wp), external :: determinant3

! Get the submatrix as four blocks, depending on whether the indices are
! greater or smaller than the element position
do ii=1,i-1
  do jj=1,j-1
    A(ii,jj) = M(ii,jj)
  end do
end do
do ii=1,i-1
  do jj=j+1,4
    A(ii,jj-1) = M(ii,jj)
  end do
end do
do ii=i+1,4
  do jj=1,j-1
    A(ii-1,jj) = M(ii,jj)
  end do
end do
do ii=i+1,4
  do jj=j+1,4
    A(ii-1,jj-1) = M(ii,jj)
  end do
end do
! The cofactor is +1/-1 times the determinant of the submatrix
f = (-One)**(i+j)
cofactor = f*determinant3(A)

end function cofactor

!***********************************************************************
!  determinant3
!
!> @brief
!>   Return the determinant of a 3&times;3 matrix.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Returns the determinant of a 3&times;3 matrix.
!>
!> @param[in] M The matrix
!>
!> @return The determinant of \p M
!***********************************************************************
function determinant3(M)

use Definitions, only: wp

implicit none
real(kind=wp) :: determinant3
real(kind=wp), intent(in) :: M(3,3)

determinant3 = (M(1,1)*M(2,2)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+M(1,3)*M(2,1)*M(3,2))- &
               (M(1,1)*M(2,3)*M(3,2)+M(1,2)*M(2,1)*M(3,3)+M(1,3)*M(2,2)*M(3,1))

end function determinant3

!***********************************************************************
!  apply_rotation
!
!> @brief
!>   Rotate a structure with a quaternion.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Apply a rotation, given by a unit quaternion, to all atoms in a structure.
!>
!> @param[in,out] x   Cartesian coordinates of the structure
!> @param[in]     nAt Number of atoms in the structure
!> @param[in]     q   A unit quaternion representing a rotation
!***********************************************************************
subroutine apply_rotation(x,nAt,q)

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: x(3,nAt)
real(kind=wp), intent(in) :: q(4)
integer(kind=iwp) :: i, iAt
real(kind=wp) :: MRot(3,3), v(3)
real(kind=wp), external :: DDot_

! Compute the rotation matrix from the quaternion
! (transposed or not, that depends on the rotation direction
!  and storage order... this version seems to work)
MRot(1,1) = q(1)**2+q(2)**2-q(3)**2-q(4)**2
MRot(2,2) = q(1)**2-q(2)**2+q(3)**2-q(4)**2
MRot(3,3) = q(1)**2-q(2)**2-q(3)**2+q(4)**2
MRot(2,1) = Two*(q(2)*q(3)+q(1)*q(4))
MRot(1,2) = Two*(q(2)*q(3)-q(1)*q(4))
MRot(3,2) = Two*(q(3)*q(4)+q(1)*q(2))
MRot(2,3) = Two*(q(3)*q(4)-q(1)*q(2))
MRot(1,3) = Two*(q(2)*q(4)+q(1)*q(3))
MRot(3,1) = Two*(q(2)*q(4)-q(1)*q(3))
! Apply the matrix to every atom in the structure
do iAt=1,nAt
  v(:) = x(:,iAt)
  do i=1,3
    x(i,iAt) = DDot_(3,MRot(:,i),1,v,1)
  end do
end do

end subroutine apply_rotation

!***********************************************************************
!> @}
!***********************************************************************
