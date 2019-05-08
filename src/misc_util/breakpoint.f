************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************

!>  @brief
!>    Set breakpoints in code.
!>
!>  @author Oskar Weser
!>
!>  @details
!>    Either you call explicitly
!>    \code{.unparsed}
!>      break breakpoint
!>    \endcode
!>    when using gdb or you add the following commands
!>    \code{.unparsed}
!>      set breakpoint pending on
!>      break breakpoint
!>    \endcode
!>    to your ~/.gdbinit.
      subroutine breakpoint()
        return
      end subroutine
