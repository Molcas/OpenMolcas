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
      Subroutine VWN_V_emb(mGrid,Rho,nRho,P2_ontop,
     &                     nP2_ontop,iSpin,F_xc,
     &                     dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,
     &                     T_X)
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*
      COMMON  / OFembed_R2/ dFMD
*
*---- Vosko-Wilk-Nusair correlation functional V
*
      Coeff=dFMD
      Call VWN_V(mGrid,Rho,nRho,iSpin,F_xc,dF_dRho,
     &           ndF_dRho,Coeff,T_X)

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
