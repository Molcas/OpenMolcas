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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2016,2017, Roland Lindh                                *
************************************************************************
      SubRoutine GMFree()
************************************************************************
*                                                                      *
*     purpose: Deallocate work space at the end of calculation to check*
*              possible errors                                         *
*                                                                      *
*     called from: SCF                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use SCF_Arrays
      use Orb_Type
      Implicit Real*8 (a-h,o-z)
*
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
#ifdef _FDE_
      ! Thomas Dresselhaus
#include "embpotdata.fh"
#endif
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      nD = iUHF + 1
*
*---- Deallocate memory
      If (Allocated(Darwin)) Then
         Call mma_deallocate(Darwin)
         Call mma_deallocate(MssVlc)
      End If
      Call mma_deallocate(KntE)
      Call mma_deallocate(EDFT)
      Call mma_deallocate(TwoHam)
      Call mma_deallocate(Vxc)
      Call mma_deallocate(Dens)
      Call mma_deallocate(OrbType)
      Call mma_deallocate(EOrb)
      Call mma_deallocate(OccNo)
      Call mma_deallocate(Fock)
      Call mma_deallocate(CMO)
      Call mma_deallocate(TrM)
*
      Call mma_deallocate(Lowdin)
      Call mma_deallocate(Ovrlp)
      Call mma_deallocate(OneHam)
      Call mma_deallocate(HDiag)
#ifdef _FDE_
      If (Allocated(Emb)) Call mma_deallocate(Emb)
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
