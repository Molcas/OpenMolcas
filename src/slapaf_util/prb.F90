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

subroutine PrB(dB,nq,nDim)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nq, nDim
real(kind=wp) :: dB(nq,nDim,nDim)
integer(kind=iwp) :: i_Dim, iq

do iq=1,nQ
  write(u6,*) ' iq=',iq
  do i_Dim=1,nDim
    write(u6,'(9F10.6)') dB(iq,i_Dim,:)
  end do
end do

return

end subroutine PrB
