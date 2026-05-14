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

module fx

! This module contains just an abstract interface, to avoid explicit interfaces

use Definitions, only: wp

implicit none
private

interface
  function f_interface(x)
    import :: wp
    real(kind=wp) :: f_interface
    real(kind=wp), intent(in) :: x
  end function f_interface
end interface

public :: f_interface

end module fx
