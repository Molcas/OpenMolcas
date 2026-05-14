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

subroutine NRed(ArrIn,ArrOut,nX,nDim,Smmtrc)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nX, nDim
real(kind=wp), intent(in) :: ArrIn(nX)
real(kind=wp), intent(out) :: ArrOut(nDim)
logical(kind=iwp), intent(in) :: Smmtrc(nX)
integer(kind=iwp) :: i_Dim, iX

i_Dim = 0
do iX=1,nX
  if (Smmtrc(iX)) then
    i_Dim = i_Dim+1
    ArrOut(i_Dim) = ArrIn(iX)
  end if
end do
if (i_Dim /= nDim) then
  write(u6,*) 'In NRed: i_Dim /= nDim'
  call Abend()
end if
!call RecPrt('ArrIn',' ',ArrIn,nX,1)
!call RecPrt('ArrOut',' ',ArrOut,nDim,1)

return

end subroutine NRed
