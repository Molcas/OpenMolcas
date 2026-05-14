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

subroutine Ext_W4hlp2(V2,M1,nc,dimab,dimapp,dimbpp,addapp,addbpp)
! this routine does:
! Extract M1(m,a",b") <- V2(m,a'b')

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nc, dimab, dimapp, dimbpp, addapp, addbpp
real(kind=wp), intent(in) :: V2(nc,dimab)
real(kind=wp), intent(out) :: M1(nc,dimapp,dimbpp)
integer(kind=iwp) :: a, ab, app

do app=1,dimapp
  a = addapp+app
  ab = nTri_Elem(a-1)+addbpp

  M1(:,app,1:dimbpp) = V2(:,ab+1:ab+dimbpp)

end do

return

end subroutine Ext_W4hlp2
