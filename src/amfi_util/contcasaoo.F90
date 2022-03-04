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
!bs contraction for powers (+2)  with alpha1*alpha3
!bs other-orbit term
!bs use averaged integrals by interchanging kinematic factors
!bs this is case a in the documentation

implicit real*8(a-h,o-z)
#include "para.fh"
#include "param.fh"
dimension ncont(4), nprim(4), primints(*), scratch1(*), scratch2(*), cont4OO(*)

ncont(1) = ncontrac(l1)
ncont(2) = ncontrac(l2)
ncont(3) = ncontrac(l3)
ncont(4) = ncontrac(l4)
nprim(1) = nprimit(l1)
nprim(2) = nprimit(l2)
nprim(3) = nprimit(l3)
nprim(4) = nprimit(l4)
ilength = nprim(1)*nprim(2)*nprim(3)*nprim(4)
nprod = ncont(1)*ncont(2)*ncont(3)*ncont(4)

!bs copy primitive integrals to scratch1
do IRUN=1,ilength
  scratch1(IRUN) = primints(IRUN)
end do
!contract   : A *alpha
!contrarray : A/E+m
!contrarray : A/E+m *alpha
!contrarray : A
!ncont      : i-th element is number of contracted functions i. index
!nprim      : i-th element is number of primitive functions  i. index
call contract(contrarray(iaddtyp2(l1)),contrarray(iaddtyp3(l2)),contrarray(iaddtyp4(l3)),contrarray(iaddtyp1(l4)),ncont,nprim, &
              scratch1,scratch2)
do irun=1,nprod
  cont4OO(nstart+irun-1) = 0.25d0*scratch1(irun)
end do

!bs copy primitive integrals to scratch1
do IRUN=1,ilength
  scratch1(IRUN) = primints(IRUN)
end do
!ncont : i-th element is number of contracted functions i. index
!nprim : i-th element is number of primitive functions  i. index
call contract(contrarray(iaddtyp4(l1)),contrarray(iaddtyp3(l2)),contrarray(iaddtyp2(l3)),contrarray(iaddtyp1(l4)),ncont,nprim, &
              scratch1,scratch2)
do irun=1,nprod
  cont4OO(nstart+irun-1) = cont4OO(nstart+irun-1)+0.25d0*scratch1(irun)
end do

!bs copy primitive integrals to scratch1
do IRUN=1,ilength
  scratch1(IRUN) = primints(IRUN)
end do
!ncont : i-th element is number of contracted functions i. index
!nprim : i-th element is number of primitive functions  i. index
call contract(contrarray(iaddtyp2(l1)),contrarray(iaddtyp1(l2)),contrarray(iaddtyp4(l3)),contrarray(iaddtyp3(l4)),ncont,nprim, &
     scratch1,scratch2)
do irun=1,nprod
  cont4OO(nstart+irun-1) = cont4OO(nstart+irun-1)+0.25d0*scratch1(irun)
end do

!bs copy primitive integrals to scratch1
do IRUN=1,ilength
  scratch1(IRUN) = primints(IRUN)
end do
!ncont : i-th element is number of contracted functions i. index
!nprim : i-th element is number of primitive functions  i. index
call contract(contrarray(iaddtyp4(l1)),contrarray(iaddtyp1(l2)),contrarray(iaddtyp2(l3)),contrarray(iaddtyp3(l4)),ncont,nprim, &
              scratch1,scratch2)
do irun=1,nprod
  cont4OO(nstart+irun-1) = cont4OO(nstart+irun-1)+0.25d0*scratch1(irun)
end do

return

end subroutine contcasaOO
