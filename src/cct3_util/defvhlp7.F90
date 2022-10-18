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

subroutine defvhlp7(r1,v,dimr1a,dimr1b,dimr1bc,dimva,dimvb,dimvc,adda)
! this routine does
! V(a,b,c)abb = R1(a,bc)
! for symb==symc
!
! r1      - r1 matrix (I)
! v       - v matrix (O)
! dimr1a  - dimension of a in R1 (I)
! dimr1b  - dimension of b (c) in R1 (I)
! dimr1bc - dimension of bc in R1 (I)
! dimva   - dimension of a in V (I)
! dimvb   - dimension of b in V (I)
! dimvc   - dimension of c in V (I)
! adda    - additional constant to a (I)

#include "t31.fh"
integer dimr1a, dimr1b, dimr1bc
integer dimva, dimvb, dimvc, adda
real*8 r1(1:dimr1a,1:dimr1bc)
real*8 v(1:dimva,1:dimvb,1:dimvc)
! help variables
integer a, b, c, bcr1

do c=1,dimvc
  do b=1,dimvb
    !bcr1 = indab(b,c)
    if (b >= c) then
      bcr1 = b*(b-1)/2+c
    else
      bcr1 = c*(c-1)/2+b
    end if
    do a=1,dimva
      v(a,b,c) = r1(a+adda,bcr1)
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(dimr1b)

end subroutine defvhlp7
