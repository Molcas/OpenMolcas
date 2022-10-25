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

subroutine defvhlp7(r1,v,dimr1a,dimr1bc,dimva,dimvb,dimvc,adda)
! this routine does
! V(a,b,c)abb = R1(a,bc)
! for symb==symc
!
! r1      - r1 matrix (I)
! v       - v matrix (O)
! dimr1a  - dimension of a in R1 (I)
! dimr1bc - dimension of bc in R1 (I)
! dimva   - dimension of a in V (I)
! dimvb   - dimension of b in V (I)
! dimvc   - dimension of c in V (I)
! adda    - additional constant to a (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr1a, dimr1bc, dimva, dimvb, dimvc, adda
real(kind=wp), intent(in) :: r1(dimr1a,dimr1bc)
real(kind=wp), intent(out) :: v(dimva,dimvb,dimvc)
integer(kind=iwp) :: b, bcr1, c

do c=1,dimvc
  do b=1,dimvb
    !bcr1 = indab(b,c)
    if (b >= c) then
      bcr1 = b*(b-1)/2+c
    else
      bcr1 = c*(c-1)/2+b
    end if
    v(:,b,c) = r1(adda+1:adda+dimva,bcr1)
  end do
end do

return

end subroutine defvhlp7
