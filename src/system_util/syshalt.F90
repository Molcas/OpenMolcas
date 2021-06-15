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
! Copyright (C) 2001, Valera Veryazov                                  *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!       replace SYSDB                                                  *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     V.Veryazov University of Lund, 2001                              *
!                                                                      *
!***********************************************************************
!  SysHalt
!
!> @brief
!>   Quit calculation
!> @author V. Veryazov
!>
!> @details
!> A routine to stop calculation because of internal error or
!> inconsistency in the code.
!>
!> @param[in] Location routine name
!***********************************************************************

subroutine SysHalt(Location)

implicit none
character(len=*), intent(in) :: Location

call SysAbendMsg(Location,'Internal error',' ')

return

end subroutine SysHalt
