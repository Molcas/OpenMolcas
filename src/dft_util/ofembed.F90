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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module OFembed

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

logical(kind=iwp) :: Do_Core = .false., Do_OFemb = .false., KEonly = .false., OFE_first = .true.
character(len=80) :: OFE_KSDFT = ''
real(kind=wp) :: dFMD = Zero, Energy_NAD, Func_A, Func_AB, Func_B, Rep_EN, ThrFThaw = Zero, V_emb, V_Nuc_AB, V_Nuc_BA, &
                 Xsigma = 1.0e4_wp
real(kind=wp), allocatable :: FMaux(:), NDSD(:,:)

public :: dFMD, Do_Core, Do_OFemb, Energy_NAD, FMaux, Func_A, Func_AB, Func_B, KEonly, NDSD, OFE_first, OFE_KSDFT, Rep_EN, &
          ThrFThaw, V_emb, V_Nuc_AB, V_Nuc_BA, Xsigma

end module OFembed
