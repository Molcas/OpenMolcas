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

subroutine MkT_C78od(T2,Tp,Tm,dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no)
! this routine does:
! T2n(be',ga',u,v) <<-
! C7                    T2+(be',ga',uv)
! C8                    T2-(be',ga',uv)
! for beSGrp/=gaSGrp

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, dimga, dimbepp, dimgapp, addbepp, addgapp, no
real(kind=wp), intent(inout) :: T2(dimbe,dimga,no,no)
real(kind=wp), intent(in) :: Tp(dimbepp,dimgapp,nTri_Elem(no)), Tm(dimbepp,dimgapp,nTri_Elem(no-1))
integer(kind=iwp) :: u, uup, uvm, uvp

! diagonal part - contribution from T+ only
do u=1,no
  uup = nTri_Elem(u)

  ! cycle over be",ga"
  T2(addbepp+1:addbepp+dimbepp,addgapp+1:addgapp+dimgapp,u,u) = T2(addbepp+1:addbepp+dimbepp,addgapp+1:addgapp+dimgapp,u,u)+ &
                                                                Tp(:,:,uup)

end do

! off diagonal - both T+,T-
uvm = 0
do u=2,no
  uvp = nTri_Elem(u-1)

  ! cycle over be",ga"
  T2(addbepp+1:addbepp+dimbepp,addgapp+1:addgapp+dimgapp,u,1:u-1) = &
    T2(addbepp+1:addbepp+dimbepp,addgapp+1:addgapp+dimgapp,u,1:u-1)+Tp(:,:,uvp+1:uvp+u-1)+Tm(:,:,uvm+1:uvm+u-1)
  T2(addbepp+1:addbepp+dimbepp,addgapp+1:addgapp+dimgapp,1:u-1,u) = &
    T2(addbepp+1:addbepp+dimbepp,addgapp+1:addgapp+dimgapp,1:u-1,u)+Tp(:,:,uvp+1:uvp+u-1)-Tm(:,:,uvm+1:uvm+u-1)

  uvm = uvm+u-1
end do

return

end subroutine MkT_C78od
