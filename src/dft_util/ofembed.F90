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

private
public :: Do_OFemb, KEonly, OFE_first, OFE_KSDFT, ThrFThaw, Xsigma, dFMD, FMaux
public :: Rep_EN, Func_AB, Func_A, Func_B, Energy_NAD, V_Nuc_AB, V_Nuc_BA, V_emb
public :: Do_Core, NDSD

logical :: Do_OFemb = .false., KEonly = .false., OFE_first = .true.
logical :: Do_Core = .false.
character(len=80) :: OFE_KSDFT = ''
real*8, allocatable :: NDSD(:,:)
real*8 :: ThrFThaw = 0.0d0, Xsigma = 1.0d4, dFMD = 0.0d0
real*8, allocatable :: FMaux(:)
real*8 :: Rep_EN, Func_AB, Func_A, Func_B, Energy_NAD, V_Nuc_AB, V_Nuc_BA, V_emb

end module OFembed
