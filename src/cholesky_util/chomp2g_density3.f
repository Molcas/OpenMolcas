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
* Copyright (C) 2010, Jonas Bostrom                                    *
************************************************************************

      SubRoutine ChoMP2g_Density3(irc,CMO)
*     Jonas Bostrom, March 2010.
*
*     Purpose: Finalize MP2 Density.

      use ChoMP2, only: MP2D, MP2W, MP2W_e, MP2D_e
      Implicit Real*8 (a-h,o-z)
      Integer irc
      Real*8 CMO(*)
#include "real.fh"
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "chomp2g.fh"

      Character(LEN=8), Parameter:: ThisNm = 'Density3'
      Character(LEN=16), Parameter:: SecNam = 'ChoMP2g_Density3'
      Integer nOccAll(8), nOrbAll(8)

      Real*8, Allocatable:: AOTriDens(:), WAOTriDens(:)

      irc=0

      Do iSym = 1, 8
         nOccAll(iSym) = nOcc(iSym) + nFro(iSym)
         nOrbAll(iSym) = nOrb(iSym) + nDel(iSym)
      End Do

      lTriDens = 0
      Do iSym = 1, nSym
         lTriDens = lTriDens + nOrbAll(iSym)*(nOrbAll(iSym)+1)/2
      End Do

      Do iSym = 1, nSym
         Do i = 1, nOrbAll(iSym)
            Do j = 1, nOrbAll(iSym)
               If((i.le. nOrb(iSym)) .and. (j.le. nOrb(iSym))) Then
                  MP2D_e(iSym)%A(i,j) = MP2D(iSym)%A(i,j)
                  MP2W_e(iSym)%A(i,j) = MP2W(iSym)%A(i,j)
               Else
                  MP2D_e(iSym)%A(i,j) = Zero
                  MP2W_e(iSym)%A(i,j) = Zero
               End If
            End Do
         End Do
      End Do
*
      Call mma_allocate( AOTriDens,lTriDens,Label=' AOTriDens')
      Call mma_allocate(WAOTriDens,lTriDens,Label='WAOTriDens')
       AOTriDens(:)=Zero
      WAOTriDens(:)=Zero
*

      Call Build_Mp2Dens( AOTriDens,lTriDens, MP2D_e,CMO,nSym,
     &                    nOrbAll, nOccAll,.True.)
      Call Build_Mp2Dens(WAOTriDens,lTriDens, MP2W_e,CMO,nSym,
     &                    nOrbAll, nOccAll,.False.)

      Call Put_D1ao_Var(AOTriDens,lTriDens)
      Call Put_Fock_Occ(WAOTriDens,lTriDens)

      Call mma_deallocate( AOTriDens)
      Call mma_deallocate(WAOTriDens)

      End
