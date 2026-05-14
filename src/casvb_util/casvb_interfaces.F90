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

module casvb_interfaces

implicit none
private

abstract interface
  subroutine ddasonc_sub( &
#                        define _CALLING_
#                        include "ddasonc_interface.fh"
                        )
    use Definitions, only: wp, iwp
#   include "ddasonc_interface.fh"
  end subroutine ddasonc_sub

  subroutine ddsol_sub( &
#                      define _CALLING_
#                      include "ddsol_interface.fh"
                      )
    use Definitions, only: wp, iwp
#   include "ddsol_interface.fh"
  end subroutine ddsol_sub

  subroutine ddres_sub( &
#                      define _CALLING_
#                      include "ddres_interface.fh"
                      )
    use Definitions, only: wp, iwp
#   include "ddres_interface.fh"
  end subroutine ddres_sub

  subroutine ddres2upd_sub( &
#                          define _CALLING_
#                          include "ddres2upd_interface.fh"
                          )
    use Definitions, only: wp, iwp
#   include "ddres2upd_interface.fh"
  end subroutine ddres2upd_sub

  subroutine opta_sub( &
#                     define _CALLING_
#                     include "opta_interface.fh"
                     )
    use Definitions, only: iwp
#   include "opta_interface.fh"
  end subroutine opta_sub

  subroutine optb_sub( &
#                     define _CALLING_
#                     include "optb_interface.fh"
                     )
    use Definitions, only: wp, iwp
#   include "optb_interface.fh"
  end subroutine optb_sub
end interface

public :: ddasonc_sub, ddsol_sub, ddres_sub, ddres2upd_sub, opta_sub, optb_sub

end module casvb_interfaces
