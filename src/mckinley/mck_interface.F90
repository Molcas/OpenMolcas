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
    use Index_Functions, only: nTri_Elem1
    use Definitions, only: wp, iwp
#   include "grd_mck_interface.fh"
  end subroutine grd_mck_kernel

  subroutine hss_kernel( &
#                       define _CALLING_
#                       include "hss_interface.fh"
                       )
    use Index_Functions, only: nTri_Elem1
    use Definitions, only: wp, iwp
#   include "hss_interface.fh"
  end subroutine hss_kernel

  subroutine oneel_mck_kernel( &
#                             define _CALLING_
#                             include "1el_mck_interface.fh"
                             )
    use Definitions, only: wp, iwp
#   include "1el_mck_interface.fh"
  end subroutine oneel_mck_kernel

  subroutine oneeldot_mck_kernel( &
#                                define _CALLING_
#                                include "1eldot_mck_interface.fh"
                                )
    use Index_Functions, only: nTri_Elem1
    use Definitions, only: wp, iwp
#   include "1eldot_mck_interface.fh"
  end subroutine oneeldot_mck_kernel

  subroutine mck_mem( &
#                    define _CALLING_
#                    include "mem_interface.fh"
                    )
    use Definitions, only: iwp
#   include "mem_interface.fh"
  end subroutine mck_mem
end interface

public :: grd_mck_kernel, hss_kernel, mck_mem, oneel_mck_kernel, oneeldot_mck_kernel

end module mck_interface
