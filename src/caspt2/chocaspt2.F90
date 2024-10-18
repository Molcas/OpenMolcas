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
Module ChoCASPT2
! Comment: numcho_pt2, InfVec_N2_PT2, and MaxVec_PT2 are copies
! of corresponding values in module Cholesky.
! Values are transferred at the beginning of caspt2.
Integer Lsplit(8),nIsplit(8),nAsplit(8),      &
        nksh(8),nkes(8),npsh(8),npes(8),          &
        numcho_pt2(8),iALGO,InfVec_N2_PT2,MaxVec_PT2,             &
        IF_CHO,NCHSPC,NHTSPC,NFTSPC,NFTSPC_TOT,                   &
        MXNVC,MXCHARR

Type ChoType
   Integer, Allocatable:: Unit(:)
   Integer, Allocatable:: ip(:)
   Integer, Allocatable:: np(:)
   Integer, Allocatable:: sp(:)
End Type ChoType

Type (ChoType) Stuff(8)

END Module ChoCASPT2
