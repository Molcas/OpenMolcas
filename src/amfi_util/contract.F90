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

subroutine contract(coeffs1,coeffs2,coeffs3,coeffs4,ncont,nprim,arr1,arr2)
!coeffs1 :  (nprim(1),ncont(1)) modified contraction coefficients
!coeffs2 :  (nprim(2),ncont(2)) modified contraction coefficients
!coeffs3 :  (nprim(3),ncont(3)) modified contraction coefficients
!coeffs4 :  (nprim(4),ncont(4)) modified contraction coefficients
!ncont   :  i-th element is number of contracted functions i. index
!nprim   :  i-th element is number of primitive functions  i. index
!arr1    :  array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
!arr2    :  array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
!bs arr1 contains at the beginning the uncontracted integrals

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: coeffs1(*), coeffs2(*), coeffs3(*), coeffs4(*)
integer(kind=iwp), intent(in) :: ncont(4), nprim(4)
real(kind=wp), intent(inout) :: arr1(*)
real(kind=wp), intent(_OUT_) :: arr2(*)
integer(kind=iwp) :: ifirst, ifourth, isec, ithird, nnew(4), nolds(4)
real(kind=wp) :: ratio1, ratio2, ratio3, ratio4, xmax

!bs makes four index transformations in a row....
!bs try to find out, which indices should be transformed first...

ratio1 = real(nprim(1),kind=wp)/real(ncont(1),kind=wp)
ratio2 = real(nprim(2),kind=wp)/real(ncont(2),kind=wp)
ratio3 = real(nprim(3),kind=wp)/real(ncont(3),kind=wp)
ratio4 = real(nprim(4),kind=wp)/real(ncont(4),kind=wp)
nolds(:) = nprim(:)
nnew(:) = nprim(:)
!bs determine first, second,third and last index
!***********************************************************************
!bs determine the first
xmax = max(ratio1,ratio2,ratio3,ratio4)
if (xmax == ratio1) then
  ifirst = 1
  ratio1 = Zero
  nnew(ifirst) = ncont(ifirst)
  call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else if (xmax == ratio2) then
  ifirst = 2
  ratio2 = Zero
  nnew(ifirst) = ncont(ifirst)
  call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else if (xmax == ratio3) then
  ifirst = 3
  ratio3 = Zero
  nnew(ifirst) = ncont(ifirst)
  call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else if (xmax == ratio4) then
  ifirst = 4
  ratio4 = Zero
  nnew(ifirst) = ncont(ifirst)
  call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else
  ifirst = 0
  write(u6,*) 'Contract: you should not be here!'
  call abend()
end if
nolds(ifirst) = nnew(ifirst)
!***********************************************************************
!bs determine the second
xmax = max(ratio1,ratio2,ratio3,ratio4)
if (xmax == ratio1) then
  isec = 1
  ratio1 = Zero
  nnew(isec) = ncont(isec)
  call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else if (xmax == ratio2) then
  isec = 2
  ratio2 = Zero
  nnew(isec) = ncont(isec)
  call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else if (xmax == ratio3) then
  isec = 3
  ratio3 = Zero
  nnew(isec) = ncont(isec)
  call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else if (xmax == ratio4) then
  isec = 4
  ratio4 = Zero
  nnew(isec) = ncont(isec)
  call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else
  isec = 0
  write(u6,*) 'Contract: you should not be here!'
  call abend()
end if
nolds(isec) = nnew(isec)
!***********************************************************************
!bs determine the third
xmax = max(ratio1,ratio2,ratio3,ratio4)
if (xmax == ratio1) then
  ithird = 1
  ratio1 = Zero
  nnew(ithird) = ncont(ithird)
  call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else if (xmax == ratio2) then
  ithird = 2
  ratio2 = Zero
  nnew(ithird) = ncont(ithird)
  call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else if (xmax == ratio3) then
  ithird = 3
  ratio3 = Zero
  nnew(ithird) = ncont(ithird)
  call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else if (xmax == ratio4) then
  ithird = 4
  ratio4 = Zero
  nnew(ithird) = ncont(ithird)
  call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
else
  ithird = 0
  write(u6,*) 'Contract: you should not be here!'
  call abend()
end if
nolds(ithird) = nnew(ithird)
!***********************************************************************
!bs determine the last
xmax = max(ratio1,ratio2,ratio3,ratio4)
if (xmax == ratio1) then
  ifourth = 1
  ratio1 = Zero
  nnew(ifourth) = ncont(ifourth)
  call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else if (xmax == ratio2) then
  ifourth = 2
  ratio2 = Zero
  nnew(ifourth) = ncont(ifourth)
  call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else if (xmax == ratio3) then
  ifourth = 3
  ratio3 = Zero
  nnew(ifourth) = ncont(ifourth)
  call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else if (xmax == ratio4) then
  ifourth = 4
  ratio4 = Zero
  nnew(ifourth) = ncont(ifourth)
  call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
else
  ifourth = 0
  write(u6,*) 'Contract: you should not be here!'
  call abend()
end if
!bs contracted integrals are now on
!bs arr1(ncont1,ncont2,ncont3,ncont4)

return

end subroutine contract
