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

module grd_interface

implicit none
private

abstract interface
  subroutine grd_kernel( &
#                       define _CALLING_
#                       include "grd_interface.fh"
                       )
    use Index_Functions, only: nTri_Elem1
    use Definitions, only: wp, iwp
#   include "grd_interface.fh"
  end subroutine grd_kernel

  subroutine grd_mem( &
#                    define _CALLING_
#                    include "mem_interface.fh"
                    )
    use Definitions, only: iwp
#   include "mem_interface.fh"
  end subroutine grd_mem
end interface

public :: grd_kernel, grd_mem

end module grd_interface
