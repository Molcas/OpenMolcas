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
! Copyright (C) 1990, IBM                                              *
!               1991, Roland Lindh                                     *
!***********************************************************************

function TF(mdc,iIrrep,iComp)

use Center_Info, only: dc
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: TF
integer(kind=iwp), intent(in) :: mdc, iIrrep, iComp
logical(kind=iwp), external :: TstFnc

TF = TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab)

end function TF
