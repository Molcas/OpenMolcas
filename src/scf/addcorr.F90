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

module AddCorr

use Definitions, only: wp, iwp

implicit none
private

logical(kind=iwp) :: Do_Addc, Do_Tw
character(len=80) :: ADDC_KSDFT
real(kind=wp) :: DE_KSDFT_c

public :: ADDC_KSDFT, DE_KSDFT_c, Do_Addc, Do_Tw

end module AddCorr
