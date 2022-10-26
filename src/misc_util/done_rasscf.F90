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

subroutine Done_RASSCF(CMO,OCC,D)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute the active one-body density                              *
!                                                                      *
!     calling arguments:                                               *
!     CMO     : input, array of real*8                                 *
!               MO-coefficients                                        *
!     OCC     : input, array of real*8                                 *
!               occupation numbers                                     *
!     D       : output, array of real*8                                *
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

implicit real*8(A-H,O-Z)
dimension CMO(*), OCC(*), D(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"

iOff1 = 0
iOff2 = 0
iOff3 = 0
do iSym=1,nSym
  iBas = nBas(iSym)
  iFro = nFro(iSym)
  iIsh = nIsh(iSym)
  iAsh = nAsh(iSym)
  if (iBas /= 0) then
    do i=1,iBas
      ii = (i*i-i)/2
      do j=1,i
        Sum = 0.0d0
        do k=1,iFro+iIsh+iAsh
          Sum = Sum+OCC(iOff3+k)*CMO(iOff1+(k-1)*iBas+i)*CMO(iOff1+(k-1)*iBas+j)
        end do
        D(iOff2+ii+j) = 2.0d0*Sum
        if (j == i) D(iOff2+ii+j) = Sum
      end do
    end do
  end if
  iOff1 = iOff1+iBas*iBas
  iOff2 = iOff2+(iBas*iBas+iBas)/2
  iOff3 = iOff3+iBas
end do

return

end subroutine Done_RASSCF
