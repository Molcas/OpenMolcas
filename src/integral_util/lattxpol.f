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
      Subroutine lattXPol(Grid,nGrid,nGrid_Eff,PolEff,DipEff,
     &                    XF,nXF,nOrd_XF,nPolComp)
*
      Implicit Real*8 (a-h,o-z)
*
      Real*8 Grid(3,nGrid), PolEff(nPolComp,nGrid), DipEff(nGrid)
      Real*8 XF(*)

#include "real.fh"

*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

*     Calculate number of entries per XFIELD point
      Inc = 3
      Do iOrdOp = 0, nOrd_XF
         Inc = Inc + nElem(iOrdOp)
      End Do
      Inc = Inc + 6  !iXpolType always .gt.0

*     Insert XFIELD polarisabilities into Grid
      Do iXF=1,nXF
         nGrid_Eff=nGrid_Eff+1
         Do j=1,nPolComp
            PolEff(j,nGrid_Eff)=XF(iXF*Inc-6+j)
         EndDo
         DipEff(nGrid_Eff)=Zero
         Do j=1,3
            Grid(j,nGrid_Eff)=XF((iXF-1)*Inc+j)
         EndDo
      EndDo

      Return
      End
