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

      SubRoutine ChoMP2g_Setup(irc,EOcc,EVir)
      use ChoMP2, only: ChoMP2g_Allocated, EFrozT, EOccuT, EVirtT
      use ChoMP2, only: AdrR1, AdrR2
      use ChoMP2, only: MP2D_full, MP2D
      use ChoMP2, only: MP2W_full, MP2W
      use ChoMP2, only: MP2D_e_full, MP2D_e
      use ChoMP2, only: MP2W_e_full, MP2W_e
      use Constants
*
*     Jonas Bostrom, Feb 2010
*
*     Purpose: Do some additional setup only needed for
*              MP2-gradients or properties.
      Implicit Real*8 (a-h,o-z)
#include "chomp2g.fh"
#include "chomp2.fh"
#include "choorb.fh"
#include "cholesky.fh"
#include "stdalloc.fh"
      Integer irc
      Real*8 EOcc(*), EVir(*)
*
******************************************************
      MulD2h(i,j)=iEor(i-1,j-1) + 1
******************************************************
      nMOType = 3
      Call ChoMP2_GetInf(nOrb,nOcc,nFro,nDel,nVir)
*
*     Initialize an  offset for writing choleskyvectors to disk
      Do iProdType = 1, nMoType**2
         Do iSym = 1, nSym
            iAdrOff(iSym,iProdType) = 0
         End Do
      End Do

      nOccVirT = nOcc(1)*nVir(1)
      Do iSym = 2, nSym
         nOccVirT = nOccVirT + nOcc(iSym)*nVir(iSym)
      End Do

      Do iMOType = 1, nMOType
         Do iSym = 1, nSym
            nMO(iSym,1) = nFro(iSym)
            nMO(isym,2) = nOcc(iSym)
            nMO(iSym,3) = nVir(iSym)
         End Do
      End Do
*
      Do iMoType = 1, nMoType
         Do jMoType = 1, nMoType
            iProdType = jMOType + (iMOType-1)*nMOType
            Do iSym = 1, nSym
               nMoMo(iSym,iProdType) = 0
               Do iSymP = 1,nSym
                  iSymQ = MulD2h(iSymP,iSym)
                  iMoMo(iSymq,iSymp,iProdType) =
     &                  nMoMo(iSym,iProdType)
                  nMoMo(iSym,iProdType) = nMoMo(iSym,iProdType) +
     &                  nMO(iSymP,iMOType)*
     &                  nMO(iSymQ,jMOType)
               End Do
            End Do
         End Do
      End Do
      Do iMoType = 1, nMoType
         Do iSym = 1, nSym
            nMoAo(iSym,iMoType) = 0
            Do iSymAl = 1, nSym
               iSymP = MulD2h(iSymAl,iSym)
               iMoAo(iSymP,iSymAl,iMoType) = nMoAo(iSym,iMoType)
               nMoAo(iSym,iMoType) =
     &               nMoAo(iSym,iMoType) +
     &               nBas(iSymAl)*nMO(iSymP,iMoType)
            End Do
         End Do
      End Do
      Do iMoType = 1, nMoType
         Do iSym = 1, nSym
            nAoMo(iSym,iMoType) = 0
            Do iSymQ = 1, nSym
               iSymAl = MulD2h(iSymQ,iSym)
               iAoMo(iSymAl,iSymQ,iMoType) = nAoMo(iSym,iMoType)
               nAoMo(iSym,iMoType) =
     &               nAoMo(iSym,iMoType) +
     &               nBas(iSymAl)*nMO(iSymQ,iMoType)
            End Do
         End Do
      End Do

*     Allocate MP2_density
*     --------------------

      lDens = nOrb(1)*nOrb(1)
      Do iSym = 2, nSym
         lDens = lDens + nOrb(iSym)*nOrb(iSym)
      End Do
*
      ChoMP2g_Allocated=.True.

      Call mma_allocate(MP2D_full,lDens,Label='MP2D_full')
      Call mma_allocate(MP2W_full,lDens,Label='MP2W_full')
      MP2D_full(:)=Zero
      MP2W_full(:)=Zero

      iE = 0
      Do iSym = 1, nSym
         nb=nOrb(iSym)
         iS = iE + 1
         iE = iE + nb**2
         MP2D(iSym)%A(1:nb,1:nb) => MP2D_full(iS:iE)
         MP2W(iSym)%A(1:nb,1:nb) => MP2W_full(iS:iE)
      End Do

*     Allocate extended MP2_density (with deleted orbitals)
*     -----------------------------------------------------

      lDens_e = (nOrb(1)+nDel(1))*(nOrb(1)+nDel(1))
      Do iSym = 2, nSym
         lDens_e = lDens_e + (nOrb(iSym)+nDel(iSym))
     &                     * (nOrb(iSym)+nDel(iSym))
      End Do

      Call mma_allocate(MP2D_e_full,lDens_e,Label='MP2D_e_full')
      Call mma_allocate(MP2W_e_full,lDens_e,Label='MP2W_e_full')
      MP2D_e_full(:)=Zero
      MP2W_e_full(:)=Zero

      iE = 0
      Do iSym = 1, nSym
         nb=nOrb(iSym)+nDel(iSym)
         iS = iE + 1
         iE = iE + nb**2
         MP2D_e(iSym)%A(1:nb,1:nb) => MP2D_e_full(iS:iE)
         MP2W_e(iSym)%A(1:nb,1:nb) => MP2W_e_full(iS:iE)
      End Do

*     Allocate adress-field for reordered R-vectors
*     ---------------------------------------------
      Call mma_allocate(AdrR1,nSym,nSym,nOccT,Label='AdrR1')
      Call mma_allocate(AdrR2,nSym,nSym,nVirT,Label='AdrR2')

*    Allocate a vector for the orbital energies of frozen and virtual
*     frozen molecules.
      Call mma_allocate(EFrozT,Max(1,nFroT),Label='EFrozT')
      Call mma_allocate(EOccuT,Max(1,nOccT),Label='EOccuT')
      Call mma_allocate(EVirtT,Max(1,nVirT),Label='EVirtT')
*     Fill them with the right things
      Do iSym = 1, nSym
         Do i = 1, nFro(iSym)
            EFrozT(iFro(iSym)+i) = EOcc(iFro(iSym)+nOccT +i)
         End Do
         Do i = 1, nOcc(iSym)
            EOccuT(iOcc(iSym)+i) =
     &                    EOcc(iOcc(iSym) + i)
         End Do
         Do i = 1, nVir(iSym)
            EVirtT(iVir(iSym)+i) =
     &                    EVir(iVir(iSym)+iDel(iSym)+ i)
         End Do
      End Do

      irc=0

      End
