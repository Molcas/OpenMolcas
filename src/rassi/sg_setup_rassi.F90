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
subroutine SG_setup_RASSI(nIrrep,NACTEl,MPLETT,SGS,CIS)

use sguga, only: CIStruct, EXStruct, SGStruct, SG_Init
use rassi_aux, only: Level
use definitions, only: iwp

integer(kind=iwp), intent(in):: nIrrep,NACTEL,MPLETT
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS

Call SG_Init(nIrrep,NACTEL,MPLETT,SGS,CIS,xLevel=Level)

End subroutine SG_setup_RASSI
