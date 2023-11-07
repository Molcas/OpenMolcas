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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine xflush(Lu)
!***********************************************************************
!                                                                      *
!   purpose:                                                           *
!   Dump buffers                                                       *
!                                                                      *
!   calling arguments:                                                 *
!   Lu      : integer                                                  *
!             Logical unit number                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!   written by:                                                        *
!   M.P. Fuelscher                                                     *
!   University of Lund, Sweden, 1996                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!   history: Flush intrinsic available since Fortran 2003 (Jun. 2021)  *
!                                                                      *
!***********************************************************************

#define _F2003_

use Definitions, only: iwp
#if ! defined (_F2003_) && defined (_CRAY_C90_)
use Definitions, only: u6
#endif

#ifdef _F2003_

implicit none
integer(kind=iwp), intent(in) :: Lu

flush(Lu)

#elif defined (_CRAY_C90_)
if (Lu == u6) then
  call flush(101)
else
  call flush(Lu)
end if
#elif defined (__INTEL_COMPILER) || defined (NAGFOR)
#include "macros.fh"
unused_var(Lu)
return
#elif defined (_SOLARIS_) || defined (_IRIX64_) || defined (_HP_UX_)
call flush(Lu)
#elif defined (_PRIMEPOWER_) || defined (_LINUX_)
call flush(Lu)
#endif

return

end subroutine xflush
