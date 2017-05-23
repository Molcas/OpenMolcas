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
* Copyright (C) 2008, Jonas Bostrom                                    *
************************************************************************
      SubRoutine ChoMP2_Dens_Setup(ip_CMO_Reord,CMO,EOcc,EVir)
C
C     Jonas Bostrom, May. 2008.
C
C     Purpose: Allocate memory and setup pointers and prepare to calculate
C              mp2-densities
C
#include "implicit.fh"
      Dimension CMO(*), EOcc(*), EVir(*)
#include "WrkSpc.fh"
#include "chomp2.fh"
#include "cholesky.fh"
#include "choorb.fh"
*
      Character*17 SecNam
      Parameter (SECNAM = 'ChoMP2_Dens_Setup')
*
      lCMO = 0
      Do iSym = 1, nSym
         lCMO = lCMO + nBas(iSym)*(nOrb(iSym) + nDel(iSym))
      End Do
*
      MaxVO = 0
      Do iSym = 1, nSym
         MaxVO = max(MaxVO,nFro(iSym) + nOcc(iSym))
         MaxVO = max(MaxVO,nVir(iSym) + nDel(iSym))
      End Do
      lTmpLvec = MaxVO*MaxVo
      Call GetMem('TmpLvec','Allo','Real',ip_tmpL,lTmpLvec)
*
      Call GetMem('CMO_reord','Allo','Real',ip_CMO_reord,lCMO)
      Call FZero(Work(ip_CMO_reord),lCMO)
*     Reordering of CMO, step 1 - cut occupied and virtual orbitals
*                                 into the matrix
      iSym_offset1 = 0
      iSym_offset2 = 0
      Do iSym = 1, nSym
         iOrb_offset = 0
         Do iOrb = 1, nOrb(iSym) + nDel(iSym)
            Call dCopy_(nBas(iSym),
     &                 CMO(iSym_offset1+iOrb_offset+1),1,
     &                 Work(ip_CMO_reord+iOrb-1+iSym_offset2),
     &                 nOrb(iSym)+nDel(iSym))
            iOrb_offset = iOrb_offset + nOrb(iSym) + nDel(iSym)
         End Do
         iSym_offset1 = iSym_offset1 + nBas(iSym)*
     &                 (nOrb(iSym)+nDel(iSym))
         iSym_offset2 = iSym_offset2 + nBas(iSym)*
     &                  (nOrb(iSym)+nDel(iSym))
      End Do
*
*
      l_Density = 0
      l_Mp2Lagr = 0
      l_DiaA    = 0
      Do iSym = 1, nSym
         l_Density = l_Density + (nOrb(iSym)+nDel(iSym))
     &                         * (nOrb(iSym)+nDel(iSym))
         l_Mp2Lagr = l_Mp2Lagr + (nFro(iSym) + nOcc(iSym))
     &                         * (nVir(iSym) + nDel(iSym))
         l_DiaA    = l_DiaA    + (nFro(iSym) + nOcc(iSym))
     &                         * (nVir(iSym) + nDel(iSym))
      End Do
      Call GetMem('MP2Density','Allo','Real',ip_First_Density,
     &            l_Density)
      Call GetMem('MP2WDensity','Allo','Real',ip_First_WDensity,
     &            l_Density)
      Call GetMem('MP2Lagr','Allo','Real',
     &            ip_First_Mp2Lagr, l_Mp2Lagr)
      Call GetMem('MP2DiaA','Allo','Real',
     &            ip_First_DiaA, l_DiaA)

      Call FZero(Work(ip_First_Density),    l_Density)
      Call FZero(Work(ip_First_WDensity),   l_Density)
      Call FZero(Work(ip_First_Mp2Lagr),    l_Mp2Lagr)
      Call FZero(Work(ip_First_DiaA),       l_DiaA)
*
      ip_Density(1)    = ip_First_Density
      ip_WDensity(1)   = ip_First_WDensity
      ip_Mp2Lagr(1)    = ip_First_Mp2Lagr
      ip_DiaA(1)       = ip_First_DiaA
      Do i = 1, nSym-1
         ip_Density(i+1)    = ip_Density(i)    + (nOrb(i)+nDel(i))
     &                                         * (nOrb(i)+nDel(i))
         ip_WDensity(i+i)   = ip_WDensity(i)   + (nOrb(i)+nDel(i))
     &                                         * (nOrb(i)+nDel(i))
         ip_Mp2Lagr(i+1)    = ip_Mp2Lagr(i)    + (nFro(i) +nOcc(i))
     &                                         * (nVir(i) +nDel(i))
         ip_DiaA(i+1)       = ip_DiaA(i)       + (nFro(i)+nOcc(i))
     &                                         * (nVir(i)+nDel(i))
      End Do
*
*    Allocate a vector for the orbital energies of frozen and virtual
*     frozen molecules.
      Call GetMem('EFro','Allo','Real',ip_EFroz,nFroT)
      Call GetMem('EOcc','Allo','Real',ip_EOccu,nOccT)
      Call GetMem('EVir','Allo','Real',ip_EVirt,nVirT)
      Call GetMem('EDel','Allo','Real',ip_EDele,nDelT)
*     Fill them with the right things
      Do iSym = 1, nSym
         Do i = 0, nFro(iSym)-1
            Work(ip_EFroz + iFro(iSym)+i) =
     &                    EOcc(iFro(iSym)+iOcc(iSym)+nOcc(iSym) +i+1)
         End Do
         Do i = 0, nOcc(iSym)-1
            Work(ip_EOccu + iOcc(iSym)+i) =
     &                    EOcc(iFro(iSym)+iOcc(iSym) + i+1)
         End Do
         Do i = 0, nVir(iSym)-1
            Work(ip_EVirt + iVir(iSym)+i) =
     &                    EVir(iVir(iSym)+iDel(iSym)+ i+1)
         End Do
         Do i=0, nDel(iSym)-1
            Work(ip_EDele + iDel(iSym)+i) =
     &                    EVir(iVir(iSym) + iDel(iSym) + nVir(iSym)+i+1)
         End Do
      End Do
      Return
      End
