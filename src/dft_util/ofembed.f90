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
Module OFembed
Private
Public:: Do_OFemb, KEonly, OFE_first, OFE_KSDFT, ThrFThaw, Xsigma, dFMD, FMaux
Public::Rep_EN,Func_AB,Func_A,Func_B,Energy_NAD,V_Nuc_AB,V_Nuc_BA,V_emb

Logical::  Do_OFemb=.False., KEonly=.False., OFE_first=.True.
Character(LEN=16)::  OFE_KSDFT=''
#ifdef _NOT_USED_
Integer:: ip_NDSD=-696696, l_NDSD=0
#endif
Real*8:: ThrFThaw=0.0D0, Xsigma=1.0d4, dFMD=0.0D0
Real*8, Allocatable:: FMaux(:)
Real*8::Rep_EN,Func_AB,Func_A,Func_B,Energy_NAD,V_Nuc_AB,V_Nuc_BA,V_emb

End Module OFembed
