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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
       Subroutine NEWF1(mGrid,Rho,nRho,P2_ontop,
     &                  nP2_ontop,iSpin,F_xc,
     &                  dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*     Call Sergey's New Functional 1                                   *
*                                                                      *
      Coeff=One
      Call Do_NewFunctional_1(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                        Coeff,iSpin,F_xc,
     &                        P2_ontop,nP2_ontop,dF_dP2ontop,
     &                        ndF_dP2ontop,T_X)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
