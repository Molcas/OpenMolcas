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
! Copyright (C) 2016, Roland Lindh                                     *
!***********************************************************************
!     This module contains the most-important globally-allocated arrays*
!     in the SCF program.                                              *
!                                                                      *
!***********************************************************************
Module SCF_Arrays
   Real*8, Dimension(:),     Allocatable:: HDiag, Ovrlp, OneHam, EDFT, KntE, Darwin, MssVlc
   Real*8, Dimension(:,:),   Allocatable:: CMO, TrM, FockAO, OccNo, EOrb, CInter, CMO_ref
   Real*8, Dimension(:,:,:), Allocatable:: TrDh, TrDP, TrDD
   Real*8, Dimension(:,:),   Allocatable, Target:: FockMO
   Real*8, Dimension(:,:,:), Allocatable, Target:: TwoHam, Vxc, Dens
End Module SCF_Arrays
