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

subroutine cartoneZ(L,Lmax,onecontr,ncontrac,MxcontL,onecartZ)
!bs arranges the cartesian one-electron integrals for Z on a square matrix

use index_functions, only: iTri
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L, Lmax, ncontrac, MxcontL
real(kind=wp), intent(in) :: onecontr(MxcontL,MxcontL,-Lmax:Lmax,3)
real(kind=wp), intent(inout) :: onecartZ(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))
integer(kind=iwp) :: iaddr, Mprime

!bs - + Integrals    m || mprime     mprime=m
do Mprime=1,L
  iaddr = iTri(Mprime+L+1,-mprime+L+1)
  onecartZ(1:ncontrac,1:ncontrac,iaddr) = onecartZ(1:ncontrac,1:ncontrac,iaddr)+ &
                                          Half*(onecontr(1:ncontrac,1:ncontrac,Mprime,2)-onecontr(1:ncontrac,1:ncontrac,-Mprime,2))
end do

return

end subroutine cartoneZ
