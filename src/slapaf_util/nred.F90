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

implicit real*8(a-h,o-z)
real*8 ArrIn(nX), ArrOut(nDim)
logical Smmtrc(nX)

iDim = 0
do iX=1,nX
  if (Smmtrc(iX)) then
    iDim = iDim+1
    ArrOut(iDim) = ArrIn(iX)
  end if
end do
if (iDim /= nDim) then
  write(6,*) 'In NRed: iDim /= nDim'
  call Abend()
end if
!call RecPrt('ArrIn',' ',ArrIn,nX,1)
!call RecPrt('ArrOut',' ',ArrOut,nDim,1)

return

end subroutine NRed
