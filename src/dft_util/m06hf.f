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
* Copyright (C) 2001, Roland Lindh                                     *
************************************************************************
      Subroutine M06HF(mGrid,Rho,nRho,P2_ontop,
     &                 nP2_ontop,iSpin,F_xc,
     &                 dF_dP2ontop,ndF_dP2ontop)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. March 2001                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
      integer ijzy
*                                                                      *
************************************************************************
*                                                                      *
      ijzy=2
      CoeffA=One*CoefX
      Call XM06(mGrid,
     &          CoeffA,iSpin,F_xc,ijzy)
      CoeffA=One*CoefX
      Call XVS98(mGrid,
     &           CoeffA,iSpin,F_xc,ijzy+1)

      CoeffA=One*CoefR
      Call CM06(mGrid,
     &          CoeffA,iSpin,F_xc,ijzy)
      CoeffA=One*CoefR
      Call CVS98(mGrid,
     &           CoeffA,iSpin,F_xc,ijzy+1)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_Integer(nRho)
         Call Unused_real_array(Rho)
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
