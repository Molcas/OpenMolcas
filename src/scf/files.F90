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
!     Define files ( file names and unit numbers )                     *
!----------------------------------------------------------------------*
!     FnOne, LuOne - logical file name and unit for one-electron       *
!                    integral file                                     *
!     FnOrd, LuOrd - logical file name and unit for two-electron       *
!                    integral file                                     *
!     FnOut, LuOut - logical file name and unit for output orbitals    *
!                    file                                              *
!     FnInp, LuInp - logical file name and unit for input orbitals     *
!                    file                                              *
!     FnDns, LuDns - logical file name and unit for input density      *
!                    file                                              *
!     FnDOu, LuDOu - logical file name and unit for output density     *
!                    file                                              *
!     FnDSt, LuDSt - logical file name and unit for file with          *
!                    density differences                               *
!     FnOSt, LuOSt - logical file name and unit for file with          *
!                    dVxcdRa contributions                             *
!     FnTSt, LuTSt - logical file name and unit for file with          *
!                    two-electron contributions                        *
!     FnGrd, LuGrd - logical file name and unit for file with          *
!                    gradients used in DIIS                            *
!     FnDGd, LuDGd - logical file name and unit for file with          *
!                    gradient differences of subsequent SO steps       *
!     Fnx, Lux     - logical file name and unit for file with          *
!                    orbital rotation parameters x (for DIIS/QNR)      *
!     FnDel, LuDel - logical file name and unit for file with          *
!                    Delta for recursive calculation of                *
!     Fny, Luy     - logical file name and unit for file with          *
!                    y-vector for recursive calculation of             *
!                    Dn=-HnGn (H:inv Hessian orb rot) 2nd order update *
!----------------------------------------------------------------------*

module Files

integer :: LuOne = 10
character(len=8), parameter :: FnOne = 'ONEINT '
integer :: LuOrd = 40
character(len=8), parameter :: FnOrd = 'ORDINT '
integer :: LuOut = 20
character(len=8), parameter :: FnOut = 'SCFORB  '
integer :: LuInp = 25
integer :: LuDSt = 34
character(len=8), parameter :: FnDSt = 'DNSMAT  '
integer :: LuOSt = 40
character(len=8), parameter :: FnOSt = 'DVXCDR  '
integer :: LuTSt = 35
character(len=8), parameter :: FnTSt = 'TWOHAM  '
integer :: LuGrd = 36
character(len=8), parameter :: FnGrd = 'GRADIENT'
integer :: LuDGd = 37
character(len=8), parameter :: FnDGd = 'SODGRAD '
integer :: Lux = 38
character(len=8), parameter :: Fnx = 'SOXVEC  '
integer :: LuDel = 39
character(len=8), parameter :: FnDel = 'SODELTA '
integer :: Luy = 29
character(len=8), parameter :: Fny = 'SOYVEC  '

end module Files
