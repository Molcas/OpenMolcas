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
* Copyright (C) 2020, Roland Lindh                                     *
************************************************************************
      Module k2_arrays
!
!     DeDe is an array with desymmetrized 1-particle densities used
!     in the integral direct construction of the Fock matrix. In this
!     the array contain the subblocks with their associatated base
!     pointer.
!     ipDeDe: the set of densities to be used sorted according to the
!      integral shell pairs they contribute. Pointers to the individual
!      matrices are stored in ipOffD.
!     ipD00: some of the pointers in ipOffD point to an "empty" slot
!      and this pointer is to this part of DeDe.
!     ipDijS: points to an auxiliary pice of memory which is used in
!      case that used a subset of the elements of a matrix is used. In
!      this case picky_ will extract those elements and put them into
!      this part of DeDe on the fly.
!
      Integer, Dimension(:,:), Allocatable :: ipOffD
      Real*8, Allocatable:: FT(:), DeDe(:)
      Integer  ipDeDe, ipD00, ipDijS
      Integer  nDeDe, nDeDe_DFT, MaxDe, MxDij, MxFT, nFT
      Logical  DoGrad_, DoHess_
      Real*8, Target, Allocatable:: Fq(:), Dq(:)
      Real*8, Pointer:: pFq(:), pDq(:)
      Integer MemR, MemI, ipZeta, ipiZet
      Real*8, Allocatable:: Mem_DBLE(:)
      Integer, Allocatable:: Mem_INT(:)
      Real*8, Allocatable:: Aux(:), Sew_Scr(:)
      Integer, Allocatable:: iSOSym(:,:)

      End Module
