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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module mck_interface

implicit none
private

abstract interface
  subroutine grd_mck_kernel( &
#                           define _CALLING_
#                           include "grd_mck_interface.fh"
                           )
    use Definitions, only: wp, iwp
#   include "grd_mck_interface.fh"
  end subroutine grd_mck_kernel

  subroutine grd_mck_mem( &
#                        define _CALLING_
#                        include "grdmem_mck_interface.fh"
                        )
    use Definitions, only: iwp
#   include "grdmem_mck_interface.fh"
  end subroutine grd_mck_mem
end interface

public :: grd_mck_kernel, grd_mck_mem

end module mck_interface
