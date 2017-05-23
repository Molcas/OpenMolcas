************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine MP2gDens_setup()
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (a-h,o-z)
*
#include "WrkSpc.fh"
#include "mp2grad.fh"
#include "corbinf.fh"
*                                                                      *
************************************************************************
*                                                                      *
      l_Density = 0
      l_Mp2Lagr = 0
      l_DiaA    = 0
      Do iSym = 1, nSym
         l_Density = l_Density + (nOrb(iSym)+nDel(iSym))
     &                         * (nOrb(iSym)+nDel(iSym))
         l_Mp2Lagr = l_Mp2Lagr + (nFro(iSym) + nOcc(iSym))
     &                         * (nExt(iSym) + nDel(iSym))
         l_DiaA    = l_DiaA    + (nFro(iSym) + nOcc(iSym))
     &                         * (nExt(iSym) + nDel(iSym))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('MP2Density','Allo','Real',ip_First_Density,
     &            l_Density)
      Call GetMem('MP2WDensity','Allo','Real',ip_First_WDensity,
     &            l_Density)
      Call GetMem('MP2Lagr','Allo','Real',
     &            ip_First_Mp2Lagr, l_Mp2Lagr)
      Call GetMem('MP2DiaA','Allo','Real',
     &            ip_First_DiaA, l_DiaA)
*
      Call FZero(Work(ip_First_Density),    l_Density)
      Call FZero(Work(ip_First_WDensity),   l_Density)
      Call FZero(Work(ip_First_Mp2Lagr),    l_Mp2Lagr)
      Call FZero(Work(ip_First_DiaA),       l_DiaA)
*                                                                      *
************************************************************************
*                                                                      *
      ip_Density(1)    = ip_First_Density
      ip_WDensity(1)   = ip_First_WDensity
      ip_Mp2Lagr(1)    = ip_First_Mp2Lagr
      ip_DiaA(1)       = ip_First_DiaA
      Do i = 1, nSym-1
         ip_Density(i+1)    = ip_Density(i)    + (nOrb(i)+nDel(i))
     &                                         * (nOrb(i)+nDel(i))
         ip_WDensity(i+1)   = ip_WDensity(i)   + (nOrb(i)+nDel(i))
     &                                         * (nOrb(i)+nDel(i))
         ip_Mp2Lagr(i+1)    = ip_Mp2Lagr(i)    + (nFro(i)+nOcc(i))
     &                                         * (nExt(i)+nDel(i))
         ip_DiaA(i+1)       = ip_DiaA(i)       + (nFro(i)+nOcc(i))
     &                                         * (nExt(i)+nDel(i))
      End Do
*
      mAdOcc(1) = ipEOcc
      nOccT = nOcc(1)
      Do iSym = 2, nSym
         mAdOcc(iSym) = mAdOcc(iSym-1) + nOcc(iSym-1)
         nOccT = nOccT + nOcc(iSym)
      End Do

      mAdVir(1) = ipEVir
      nExtT = nExt(1)
      Do iSym=2, nSym
         mAdVir(iSym) = mAdVir(iSym-1) + nExt(iSym-1)
         nExtT = nExtT + nExt(iSym)
      End Do
*
      mAdFro(1) = ipEOcc+nOccT
      Do iSym = 2, nSym
         mAdFro(iSym) = mAdFro(iSym-1) + nFro(iSym-1)
      End Do
*
      mAdDel(1) = ipEVir+nExtT
      Do iSym = 2, nSym
         mAdDel(iSym) = mAdDel(iSym-1) + nDel(iSym-1)
      End Do
*
      Return
      End
