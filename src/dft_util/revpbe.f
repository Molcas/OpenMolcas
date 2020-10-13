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
* Copyright (C) 2005, Per Ake Malmqvist                                *
************************************************************************
      Subroutine REVPBE(mGrid,Rho,nRho,P2_ontop,
     &               nP2_ontop,iSpin,F_xc,
cGLM     &               nP2_ontop,iSpin,F_xc,F_xca,F_xcb,
     &               dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,
     &               T_X)
************************************************************************
*                                                                      *
* Object: To compute the sum of the functional c_pbe and x_pbe as      *
* described in the Density Functional Repository,                      *
* (http://www.cse.clrc.ac.uk/qcg/dft) and in the article               *
*    J.P. Perdew, K. Burke, and M. Ernzerhof,                          *
*  "Generalized gradient approximation made simple,"                   *
*      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. December 2005                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
cGLM     &       F_xca(mGrid),F_xcb(mGrid),tmpB(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *

      CoeffB=One*CoefX
      Call XrevPBE(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          CoeffB,iSpin,F_xc,T_X)

      CoeffA=One*CoefR
      Call CPBE(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          CoeffA,iSpin,F_xc,T_X)
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
