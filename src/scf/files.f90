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
!
Module Files
Integer                     :: LuOne = 10
Character(LEN=8), Parameter :: FnOne = 'ONEINT '
Integer                     :: LuOrd = 40
Character(LEN=8), Parameter :: FnOrd = 'ORDINT '
Integer                     :: LuOut = 20
Character(LEN=8), Parameter :: FnOut = 'SCFORB  '
Integer                     :: LuInp = 25
Integer                     :: LuDSt = 34
Character(LEN=8), Parameter :: FnDSt = 'DNSMAT  '
Integer                     :: LuOSt = 40
Character(LEN=8), Parameter :: FnOSt = 'DVXCDR  '
Integer                     :: LuTSt = 35
Character(LEN=8), Parameter :: FnTSt = 'TWOHAM  '
Integer                     :: LuGrd = 36
Character(LEN=8), Parameter :: FnGrd = 'GRADIENT'
Integer                     :: LuDGd = 37
Character(LEN=8), Parameter :: FnDGd = 'SODGRAD '
Integer                     :: Lux   = 38
Character(LEN=8), Parameter :: Fnx   = 'SOXVEC  '
Integer                     :: LuDel = 39
Character(LEN=8), Parameter :: FnDel = 'SODELTA '
Integer                     :: Luy   = 29
Character(LEN=8), Parameter :: Fny   = 'SOYVEC  '
End Module Files
