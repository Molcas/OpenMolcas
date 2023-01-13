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
! Copyright (C) 1999, Thorstein Thorsteinsson                          *
!***********************************************************************

subroutine Real2Int(RealArg,IntArg)
!***********************************************************************
!                                                                      *
!    Purpose: Convert Real argument to integer                         *
!                                                                      *
!    Calling parameters:                                               *
!    RealArg: contains real to be converted                            *
!    IntArg : contains on return the corresponding integer array       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     T. Thorsteinsson                                                 *
!     University of Lund, Sweden, 1999                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, RtoI

implicit none
real(kind=wp), intent(in) :: RealArg
integer(kind=iwp), intent(out) :: IntArg(RtoI)

IntArg(:) = transfer(RealArg,IntArg)

return

end subroutine Real2Int
