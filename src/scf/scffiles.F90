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
!----------------------------------------------------------------------*
! Define files ( file names and unit numbers )                         *
!----------------------------------------------------------------------*
! FnOne, LuOne - logical file name and unit for one-electron integral  *
!                file                                                  *
! FnOrd, LuOrd - logical file name and unit for two-electron integral  *
!                file                                                  *
! FnOut, LuOut - logical file name and unit for output orbitals file   *
! FnInp, LuInp - logical file name and unit for input orbitals file    *
! FnDns, LuDns - logical file name and unit for input density file     *
! FnDOu, LuDOu - logical file name and unit for output density file    *
! FnDSt, LuDSt - logical file name and unit for file with density      *
!                differences                                           *
! FnOSt, LuOSt - logical file name and unit for file with dVxcdRa      *
!                contributions                                         *
! FnTSt, LuTSt - logical file name and unit for file with two-electron *
!                contributions                                         *
! FnGrd, LuGrd - logical file name and unit for file with gradients    *
!                used in DIIS                                          *
! FnlGd, LulGd - logical file name and unit for file with local grad.  *
!                used in DIIS                                          *
! FnDGd, LuDGd - logical file name and unit for file with gradient     *
!                differences of subsequent SO steps                    *
! Fnx, Lux     - logical file name and unit for file with orbital      *
!                rotation parameters x (for DIIS/QNR)                  *
! FnDel, LuDel - logical file name and unit for file with Delta for    *
!                recursive calculation of                              *
! Fny, Luy     - logical file name and unit for file with y-vector for *
!                recursive calculation of Dn=-HnGn (H:inv Hessian orb  *
!                rot) 2nd order update                                 *
!----------------------------------------------------------------------*

module SCFFiles

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: LuDel = 39, LuDGd = 37, LuDSt = 34, LuGrd = 36, LuInp = 25, LulGd = 41, LuOne = 10, LuOrd = 40, LuOSt = 40, &
                     LuOut = 20, LuTSt = 35, Lux = 38, Luy = 29
character(len=8), parameter :: FnDel = 'SODELTA', FnDGd = 'SODGRAD', FnDSt = 'DNSMAT', FnGrd = 'GRADIENT', FnlGd = 'LOCGRAD', &
                               FnOne = 'ONEINT', FnOrd = 'ORDINT', FnOSt = 'DVXCDR', FnOut = 'SCFORB', FnTSt = 'TWOHAM', &
                               Fnx = 'SOXVEC', Fny = 'SOYVEC'

public :: FnDel, FnDGd, FnDSt, FnGrd, FnlGd, FnOne, FnOrd, FnOSt, FnOut, FnTSt, Fnx, Fny, LuDel, LuDGd, LuDSt, LuGrd, LuInp, &
          LulGd, LuOne, LuOrd, LuOSt, LuOut, LuTSt, Lux, Luy

end module SCFFiles
