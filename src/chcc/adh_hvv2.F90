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

subroutine AdH_Hvv2(H,Hvv,dima,dimbe,adda,addbe,nv)
! this routine does:
! Hvv(a,be) <<- - H(be',a')

implicit none
integer dima, dimbe, nv, adda, addbe
real*8 H(1:dimbe,1:dima)
real*8 Hvv(1:nv,1:nv)
! help variables
integer a, be

do be=1,dimbe
  do a=1,dima
    Hvv(adda+a,addbe+be) = Hvv(adda+a,addbe+be)-H(be,a)
  end do
end do

return

end subroutine AdH_Hvv2
