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
*
*     Jonas Bostrom, Feb 2010
*
*     Purpose: Do some additional setup only needed for
*              MP2-gradients or properties.
#include "implicit.fh"
#include "WrkSpc.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "choorb.fh"
#include "cholesky.fh"

      Dimension EOcc(*), EVir(*)
*
******************************************************
      MulD2h(i,j)=iEor(i-1,j-1) + 1
*
      iMoMoTable(iOrb,iSym,iPar) = ipMoMoTable +
     &                             iPar-1 + (iSym-1)*3 +
     &                             (iOrb-1)*nSym*3
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


*     Stupid memoryjumping here that should be solved
      lMoMoTable = nOccVirT*3*nSym
      Call GetMem('MoMoTable','Allo','Inte',ipMoMoTable,lMoMoTable)
      Do iSym = 1, nSym
         index = 0
         Do iSymI = 1, nSym
            iSymA = MulD2h(iSymI,iSym)
            Do iI = 1, nOcc(iSymI)
               Do iA = 1, nVir(iSymA)
                  index = index + 1
                  iWork(iMoMoTable(index,iSym,1)) = iSymI
                  iWork(iMoMoTable(index,iSym,2)) = iI
                  iWork(iMoMoTable(index,iSym,3)) = iA
               End Do
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
      Call GetMem('MP2Density','Allo','Real',ipMP2D, lDens)
      Call GetMem('MP2WDensity','Allo','Real',ipMP2W, lDens)
      Call FZero(Work(ipMP2D),lDens)
      Call FZero(Work(ipMP2W),lDens)
      ipDensity(1) = ipMP2D
      ipWDensity(1) = ipMP2W
      Do i = 2, nSym
         ipDensity(i) = ipDensity(i-1) + nOrb(i-1)*nOrb(i-1)
         ipWDensity(i) = ipWDensity(i-1) + nOrb(i-1)*nOrb(i-1)
      End Do

*     Allocate extended MP2_density (with deleted orbitals)
*     -----------------------------------------------------

      lDens_e = (nOrb(1)+nDel(1))*(nOrb(1)+nDel(1))
      Do iSym = 2, nSym
         lDens_e = lDens_e + (nOrb(iSym)+nDel(iSym))
     &                     * (nOrb(iSym)+nDel(iSym))
      End Do
      Call GetMem('MP2Density_e','Allo','Real',ipMP2D_e, lDens_e)
      Call GetMem('MP2WDensity_e','Allo','Real',ipMP2W_e, lDens_e)
      Call FZero(Work(ipMP2D_e),lDens_e)
      Call FZero(Work(ipMP2W_e),lDens_e)
      ipDensity_e(1) = ipMP2D_e
      ipWDensity_e(1) = ipMP2W_e
      Do i = 2, nSym
         ipDensity_e(i) = ipDensity_e(i-1) + (nOrb(iSym)+nDel(iSym))
     &                                     * (nOrb(iSym)+nDel(iSym))
         ipWDensity_e(i) = ipWDensity_e(i-1) + (nOrb(iSym)+nDel(iSym))
     &                                       * (nOrb(iSym)+nDel(iSym))
      End Do


*     Allocate adress-field for reordered R-vectors
*     ---------------------------------------------
      lAdrR1 = nSym*nSym*nOccT
      lAdrR2 = nSym*nSym*nVirT
      Call GetMem('AdrVector1','Allo','Inte',ipAdrR1, lAdrR1)
      Call GetMem('AdrVector2','Allo','Inte',ipAdrR2, lAdrR2)

*    Allocate a vector for the orbital energies of frozen and virtual
*     frozen molecules.
      Call GetMem('EFro','Allo','Real',ip_EFroz,nFroT)
      Call GetMem('EOcc','Allo','Real',ip_EOccu,nOccT)
      Call GetMem('EVir','Allo','Real',ip_EVirt,nVirT)
*     Fill them with the right things
      Do iSym = 1, nSym
         Do i = 0, nFro(iSym)-1
            Work(ip_EFroz + iFro(iSym)+i) =
     &                    EOcc(iFro(iSym)+nOccT +i+1)
         End Do
         Do i = 0, nOcc(iSym)-1
            Work(ip_EOccu + iOcc(iSym)+i) =
     &                    EOcc(iOcc(iSym) + i+1)
         End Do
         Do i = 0, nVir(iSym)-1
            Work(ip_EVirt + iVir(iSym)+i) =
     &                    EVir(iVir(iSym)+iDel(iSym)+ i+1)
         End Do
      End Do
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(irc)
      End
