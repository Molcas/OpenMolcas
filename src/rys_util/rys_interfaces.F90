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

module Rys_Interfaces

use Definitions, only: wp, iwp

implicit none
private

abstract interface
  subroutine cff2d_kernel( &
#                         define _CALLING_
#                         include "cff2d_interface.fh"
                         )
    import :: wp, iwp
#   include "cff2d_interface.fh"
  end subroutine cff2d_kernel

  subroutine modu2_kernel( &
#                         define _CALLING_
#                         include "modu2_interface.fh"
                         )
    import :: wp, iwp
#   include "modu2_interface.fh"
  end subroutine modu2_kernel

  subroutine rys2d_kernel( &
#                         define _CALLING_
#                         include "rys2d_interface.fh"
                         )
    import :: wp, iwp
#   include "rys2d_interface.fh"
  end subroutine rys2d_kernel

  subroutine tval_kernel( &
#                        define _CALLING_
#                        include "tval_interface.fh"
                         )
    import :: wp, iwp
#   include "tval_interface.fh"
  end subroutine tval_kernel

  subroutine tval1_kernel( &
#                         define _CALLING_
#                         include "tval1_interface.fh"
                         )
    import :: wp, iwp
#   include "tval1_interface.fh"
  end subroutine tval1_kernel
end interface

public :: cff2d_kernel, modu2_kernel, rys2d_kernel, tval_kernel, tval1_kernel

end module Rys_Interfaces
