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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Tra2A(ij_pair,ij_Bas_pairs,kl_Orb_pairs,
     &                 kSym,lSym,
     &                 kBas,lBas,
     &                 kAsh,lAsh,
     &                 CMO_k,CMO_l,
     &                 IJKL,IJKX,IJVX,VXIJ)
************************************************************************
*                                                                      *
*     run the first two quarter AO --> MO transformations with         *
*     both transformed indices being active, i.e.                      *
*     (ij!kl) --> (ij!vx)                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 CMO_k(kBas,kAsh),  CMO_l(lBas,lAsh),
     &       IJKL(lBas*kBas), IJKX(kBas*lAsh),
     &       IJVX(lAsh*kAsh), VXIJ(ij_bas_pairs,kl_Orb_pairs)
*
*     Call qEnter('Tra2A')
*
*     (ij|kl) -> (ij|kx)
      Call DGEMM_('T','N',
     &            kBas,lAsh,lBas,
     &            1.0d0,IJKL,lBas,
     &            CMO_l,lBas,
     &            0.0d0,IJKX,kBas)
*
*     (ij|kx) -> (ij|vx)
      If ( kSym.eq.lSym ) then
*------ Triangular storage of target
        Call MxMt(IJKX,kBas,1,
     &            CMO_k,1,kBas,
     &            IJVX,kAsh,kBas)
      Else
*------ Rectangular storage of target
        Call DGEMM_('T','N',
     &              lAsh,kAsh,kBas,
     &              1.0d0,IJKX,kBas,
     &              CMO_k,kBas,
     &              0.0d0,IJVX,lAsh)
      End If
*
*     Scatter result onto half-transformed integral
*     list
      Do kl = 1,kl_Orb_pairs
        VXIJ(ij_pair,kl) = IJVX(kl)
      End Do
*
*     Call qExit('Tra2A')
*
      Return
      End
