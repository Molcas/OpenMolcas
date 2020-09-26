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
* Copyright (C) 2017, Giovanni Li Manni                                *
*               2017, Aron Cohen                                       *
************************************************************************
      Subroutine S12h(mGrid,Rho,nRho,P2_ontop,
     &                nP2_ontop,iSpin,F_xc,
     &                dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,
     &                T_X)
************************************************************************
*                                                                      *
* Object:     Combination of exchange S12h (25% HF exchange) and PBE   *
*             correlation terms                                        *
*             Reference:  Swart, Marcel                                *
*             Chemical Physics Letters 580 (2013) 166-171              *
*                                                                      *
*      Author: G. Li Manni and A. Cohen, Department of Electronic      *
*              Structure Theory, Max Planck Institute, Stuttgart       *
*              2017                                                    *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"

      Real*8 Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
      Integer gh_switch
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange
*
C      Coeff=1.079966d0
C
C      Call Diracx(mGrid,Rho,nRho,iSpin,F_xc,
C     &            dF_dRho,ndF_dRho,Coeff,T_X)
*
*---- S12h has its LDA part included!
*
      Coeff=0.75d0
      gh_switch = 2
      Call xS12gh(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          Coeff,iSpin,F_xc,T_X,gh_switch)
*
*---- PBE Correlation
*
      Coeff=1.0d0
      Call CPBE(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &         Coeff,iSpin,F_xc,T_X)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
