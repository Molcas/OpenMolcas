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

function Error_for_t(t,rMP,xrMP,xxrMP,xnrMP,EC,A,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Error_for_t
integer(kind=iwp), intent(in) :: ij, l, nij, lMax, nElem, nAtoms, nPert, iPrint_Errors
real(kind=wp), intent(in) :: t, rMP(nij,0:nElem-1,0:nPert-1), xnrMP(nij,nElem), EC(3,nij), R_ij(3), C_o_C(3)
real(kind=wp), intent(out) :: xrMP(nij,nElem), xxrMP(nij,nElem), Scratch_New(nij*(2+lMax+1)), Scratch_Org(nij*(2+lMax+1))
real(kind=wp), intent(inout) :: A(3,nij)
integer(kind=iwp) :: ij_temp
real(kind=wp) :: Error

A(1,ij) = EC(1,ij)+t*R_ij(1)
A(2,ij) = EC(2,ij)+t*R_ij(2)
A(3,ij) = EC(3,ij)+t*R_ij(3)
xrMP(:,:) = rMP(:,:,0)
do ij_temp=1,nij
  call ReExpand(xrMP,nij,nElem,EC(1,ij_temp),C_o_C,ij_temp,lMax)
end do
xrMP(:,:) = xrMP(:,:)+xnrMP(:,:)
xxrMP(:,:) = xrMP(:,:)
call CutOff_Error(l,lMax,xrMP,xxrMP,nij,A,C_o_C,nElem,Scratch_New,Scratch_Org,nAtoms,iPrint_Errors,Error)

Error_for_t = Error

return

end function Error_for_t
