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

function gammln(x)
! This is a Lanczos approximation for ln(Gamma(x))
! The constants below appear in multiple sources, but few give a reason
!
! The values can be computed with the following Maxima script:
! (http://www.mrob.com/pub/ries/lanczos-gamma.html)
!=======================================================================
! load("diag");
! Dc(n) := diag(makelist(2*double_factorial(2*k-1),k,0,n));
! cmatrix_element[row,col]:=
!   if is(col>row) then 0
!   else if row=1 and col=1 then 1/2
!   else (-1)^(row+col)*4^(col-1)*(row-1)*(row+col-3)!/(row-col)!/(2*col-2)!;
! C(n) := genmatrix(cmatrix_element, n+1);
! f(g,n):=sqrt(2)*(%e/(2*(n+g)+1))^(n+1/2);
! Dr(k) := diag(append([1],makelist(-(2*n+2)!/(2*n!*(n+1)!),n,0,k-1)));
! bmatrix_element[row,col] :=
!   if row = 1 then 1
!   else if is(row > col) then 0
!   else (-1)^(col-row)*binomial(col+row-3,2*row-3);
! B(k) := genmatrix(bmatrix_element,k+1);
! lanczos_coeff(g, n) :=
!   block([M : (Dr(n) . B(n)) . (C(n) . Dc(n)),
!          f : transpose(matrix(makelist(f(g,k), k, 0, n)))],
!     (M . f));
! /* Remove bfloat for symbolic expressions */
! tlg1(g,n) := bfloat(lanczos_coeff(g,n-1)*exp(g)/sqrt(2*%pi));
!
! fpprec:100;
! fpprintprec:20;
! print(tlg1(5,7));
! quit();
!=======================================================================

use Constants, only: One, Two, Half, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: gammln
real(kind=wp), intent(in) :: x
integer(kind=iwp) :: j
real(kind=wp) :: ser, tmp, y
integer(kind=iwp), parameter :: g = 5, n = 7
real(kind=wp), parameter :: cof(0:n-1) = [1.000000000190014824_wp,76.180091729471463483_wp,-86.505320329416767652_wp, &
                                          24.01409824083091049_wp,-1.2317395724501553875_wp,1.2086509738661785061e-3_wp, &
                                          -5.3952393849531283785e-6_wp], &
                            stp = sqrt(Two*Pi)

y = x
tmp = x+real(g,kind=wp)+Half
tmp = (x+Half)*log(tmp)-tmp
ser = cof(0)
do j=1,n-1
  y = y+One
  ser = ser+cof(j)/y
end do
gammln = tmp+log(stp*ser/x)

return

end function gammln
