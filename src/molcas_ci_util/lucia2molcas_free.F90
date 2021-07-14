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

subroutine LUCIA2MOLCAS_FREE()

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: LCONF, LDET
#include "csfbas.fh"

! fake values
LCONF = 1
LDET = 1

call GetMem('KICONF','Free','Integer',KICONF(1),LCONF)
call GetMem('KICTS','Free','Integer',KICTS(1),LDET)

end subroutine LUCIA2MOLCAS_FREE
