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

subroutine GATskL_Zero(igaTsk)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: igaTsk

#ifdef _MOLCAS_MPP_
call ga_zero(igaTsk)
#else
#include "macros.fh"
unused_var(igaTsk)
#endif

end subroutine GATskL_Zero
