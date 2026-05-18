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

subroutine Set_Print_Level()

use caspt2_global, only: iPrGlb
use PrintLevel, only: SILENT, USUAL
use Definitions, only: iwp

implicit none
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

! Print levels request from molcas?
IPRGLB = IPRINTLEVEL(-1)
! If inside an optimization loop, minimize output
! unless we *really* want a lot of output.
if (Reduce_Prt()) IPRGLB = max(IPRGLB-USUAL,SILENT)

end subroutine Set_Print_Level
