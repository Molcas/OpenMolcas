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

subroutine Int2Real(IntArg,RealArg)
!***********************************************************************
!                                                                      *
!    Purpose: Convert integer argument to Real                         *
!                                                                      *
!    Calling parameters:                                               *
!    IntArg : contains integer array to be converted                   *
!    RealArg: contains on return the corresponding real                *
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
integer(kind=iwp), intent(in) :: IntArg(RtoI)
real(kind=wp), intent(out) :: RealArg

RealArg = transfer(IntArg,RealArg)

return

end subroutine Int2Real
