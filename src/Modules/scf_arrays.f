************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2016, Roland Lindh                                     *
************************************************************************
*     This module contains the most-important globally-allocated arrays*
*     in the SCF program.                                              *
*                                                                      *
************************************************************************
      Module SCF_Arrays
#ifdef _FDE_
      Real*8, Dimension(:), Allocatable:: Emb
#endif
      Real*8, Dimension(:), Allocatable:: Ovrlp, OneHam, EDFT, KntE,
     &                                    Darwin, MssVlc,CorPA
      Real*8, Dimension(:,:,:), Allocatable:: TwoHam, Vxc, Dens
      Real*8, Dimension(:,:), Allocatable:: CMO, TrM, Fock, Lowdin,
     &                                      OccNo, EOrb, HDiag
      Real*8, Dimension(:,:), Allocatable:: TrDh, TrDP, TrDD, CInter
      End Module SCF_Arrays
