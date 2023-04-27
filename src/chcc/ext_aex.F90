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

subroutine Ext_Aex(Aex,VV,no)
! this routine defs:
! VV(i,u,v,j) <- Aex(ij,u,v)
! for Aex(i,j,u,v) = Aex(j,i,v,u), Aex stored only for i>=j

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: no
real(kind=wp), intent(in) :: Aex(nTri_Elem(no),no,no)
real(kind=wp), intent(out) :: VV(no,no,no,no)
integer(kind=iwp) :: i, ij, u, v

do u=1,no
  do v=1,no
    ij = 0
    do i=1,no

      VV(i,u,v,:) = Aex(ij+1:ij+i,u,v)
      VV(:,v,u,i) = Aex(ij+1:ij+i,u,v)

      ij = ij+i
    end do
  end do
end do

return

end subroutine Ext_Aex
