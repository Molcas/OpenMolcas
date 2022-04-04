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

subroutine contcasaOO(l1,l2,l3,l4,nstart,primints,scratch1,scratch2,cont4OO)
!bs contraction for powers (+2) with alpha1*alpha3
!bs other-orbit term
!bs use averaged integrals by interchanging kinematic factors
!bs this is case a in the documentation

use AMFI_global, only: contrarray, ncontrac, nprimit
use Constants, only: Quart
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: l1, l2, l3, l4, nstart
real(kind=wp), intent(in) :: primints(*)
real(kind=wp), intent(_OUT_) :: scratch1(*), scratch2(*), cont4OO(*)
integer(kind=iwp) :: ilength, ncont(4), nprim(4), nprod

ncont(1) = ncontrac(l1)
ncont(2) = ncontrac(l2)
ncont(3) = ncontrac(l3)
ncont(4) = ncontrac(l4)
nprod = ncont(1)*ncont(2)*ncont(3)*ncont(4)
nprim(1) = nprimit(l1)
nprim(2) = nprimit(l2)
nprim(3) = nprimit(l3)
nprim(4) = nprimit(l4)
ilength = nprim(1)*nprim(2)*nprim(3)*nprim(4)

!bs copy primitive integrals to scratch1
scratch1(1:ilength) = primints(1:ilength)
!contract   : A *alpha
!contrarray : A/E+m
!contrarray : A/E+m *alpha
!contrarray : A
!ncont      : i-th element is number of contracted functions i. index
!nprim      : i-th element is number of primitive functions  i. index
call contract(contrarray(:,2,l1),contrarray(:,3,l2),contrarray(:,4,l3),contrarray(:,1,l4),ncont,nprim,scratch1,scratch2)
cont4OO(nstart:nstart+nprod-1) = Quart*scratch1(1:nprod)

!bs copy primitive integrals to scratch1
scratch1(1:ilength) = primints(1:ilength)
!ncont : i-th element is number of contracted functions i. index
!nprim : i-th element is number of primitive functions  i. index
call contract(contrarray(:,4,l1),contrarray(:,3,l2),contrarray(:,2,l3),contrarray(:,1,l4),ncont,nprim,scratch1,scratch2)
cont4OO(nstart:nstart+nprod-1) = cont4OO(nstart:nstart+nprod-1)+Quart*scratch1(1:nprod)

!bs copy primitive integrals to scratch1
scratch1(1:ilength) = primints(1:ilength)
!ncont : i-th element is number of contracted functions i. index
!nprim : i-th element is number of primitive functions  i. index
call contract(contrarray(:,2,l1),contrarray(:,1,l2),contrarray(:,4,l3),contrarray(:,3,l4),ncont,nprim,scratch1,scratch2)
cont4OO(nstart:nstart+nprod-1) = cont4OO(nstart:nstart+nprod-1)+Quart*scratch1(1:nprod)

!bs copy primitive integrals to scratch1
scratch1(1:ilength) = primints(1:ilength)
!ncont : i-th element is number of contracted functions i. index
!nprim : i-th element is number of primitive functions  i. index
call contract(contrarray(:,4,l1),contrarray(:,1,l2),contrarray(:,2,l3),contrarray(:,3,l4),ncont,nprim,scratch1,scratch2)
cont4OO(nstart:nstart+nprod-1) = cont4OO(nstart:nstart+nprod-1)+Quart*scratch1(1:nprod)

return

end subroutine contcasaOO
