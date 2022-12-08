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
! Copyright (C) 1993,1999, Markus P. Fuelscher                         *
!               1993, Per Ake Malmqvist                                *
!***********************************************************************

subroutine I4Len(nData,InBuf,lIndx)
!***********************************************************************
!                                                                      *
!     purpose: return the length of packed integer values              *
!                                                                      *
!    Calling parameters:                                               *
!    nData : number of unpacked integrals                              *
!    InBuf : contains on input unpack integers                         *
!    lIndx : contains length of packed integers on output              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and P. A. Malmqvist                               *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     replace various older algorithms                                 *
!     M.P. Fuelscher, University of Lund, Sweden, 1999                 *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nData, InBuf(nData)
integer(kind=iwp), intent(out) :: lIndx(nData)
integer(kind=iwp) :: iOpt

!----------------------------------------------------------------------*

iOpt = 1
call iziplen(iOpt,nData,InBuf,lIndx)

!----------------------------------------------------------------------*

return

end subroutine I4Len
