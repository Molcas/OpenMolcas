!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Done_CASPT2(CMO,nCMO,OCC,nOCC,D,nD)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute the active one-body density                              *
!                                                                      *
!     calling arguments:                                               *
!     CMO     : input, array of real(kind=wp)                          *
!               MO-coefficients                                        *
!     OCC     : input, array of real(kind=wp)                          *
!               occupation numbers                                     *
!     D       : output, array of real(kind=wp)                         *
!               total one-body density                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use caspt2_module, only: nSym, nBas
use Constants, only: Zero, Two
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: nCMO, nOcc, nD
real(kind=wp), intent(in) :: CMO(nCMO), OCC(nOCC)
real(kind=wp), intent(out) :: D(nD)
integer(kind=iwp) iOff1, iOff2, iOff3, iSym, iBas, i, ii, j, k
real(kind=wp) :: Sum

iOff1 = 0
iOff2 = 0
iOff3 = 0
do iSym=1,nSym
  iBas = nBas(iSym)
  if (iBas == 0) cycle
  do i=1,iBas
    ii = (i*i-i)/2
    do j=1,i
      Sum = Zero
      do k=1,iBas
        Sum = Sum+OCC(iOff3+k)*CMO(iOff1+(k-1)*iBas+i)*CMO(iOff1+(k-1)*iBas+j)
      end do
      D(iOff2+ii+j) = Two*Sum
      if (j == i) D(iOff2+ii+j) = Sum
    end do
  end do
  iOff1 = iOff1+iBas*iBas
  iOff2 = iOff2+(iBas*iBas+iBas)/2
  iOff3 = iOff3+iBas
end do

end subroutine Done_CASPT2
