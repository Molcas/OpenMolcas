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

subroutine Dmat_MCLR(CMO,OCC,D)
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

use Index_Functions, only: nTri_Elem
use input_mclr, only: nBas, nSym
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*), OCC(*)
real(kind=wp), intent(_OUT_) :: D(*)
integer(kind=iwp) :: i, iBas, ii, iOff1, iOff2, iOff3, iSym, j, k
real(kind=wp) :: rSum

iOff1 = 0
iOff2 = 0
iOff3 = 0
do iSym=1,nSym
  iBas = nBas(iSym)
  if (iBas /= 0) then
    do i=1,iBas
      ii = nTri_Elem(i-1)
      do j=1,i
        rSum = Zero
        do k=1,ibas
          rSum = rSum+OCC(iOff3+k)*CMO(iOff1+(k-1)*iBas+i)*CMO(iOff1+(k-1)*iBas+j)
        end do
        if (j == i) then
          D(iOff2+ii+j) = rSum
        else
          D(iOff2+ii+j) = Two*rSum
        end if
      end do
    end do
  end if
  iOff1 = iOff1+iBas*iBas
  iOff2 = iOff2+nTri_Elem(iBas)
  iOff3 = iOff3+iBas
end do

end subroutine Dmat_MCLR
