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

subroutine Oldge(iAcc,Etot,Eold,Ract,Rold)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iAcc
real(kind=wp) :: Etot, Eold, Ract, Rold
#include "maxi.fh"
#include "qminp.fh"
integer(kind=iwp) :: i, icCom, j, k

iAcc = iAcc-1
Etot = Eold
Ract = Rold
icCom = 0
do i=1,nPart
  do j=1,nCent
    icCom = icCom+1
    do k=1,3
      Cordst(icCom,k) = OldGeo(icCom,k)
    end do
  end do
end do

return

end subroutine Oldge
