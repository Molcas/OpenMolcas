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
*               2010, Grigory A. Shamov                                *
************************************************************************
      Subroutine SSB(mGrid,Rho,nRho,P2_ontop,
     &                nP2_ontop,iSpin,F_xc,
     &                dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,
     &                T_X)
************************************************************************
*                                                                      *
* Object:     Combination of excnange SSB-sw and PBE correlation terms *
*             ref (secondary): Swart, Sola, Bickelhaupt                *
*             J.Chem.Phys. 131 (2009) 094103. Note that it isnt SSB-D  *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. March 2001                              *
*             Grigory A Shamov, U of Mantoba, 2010                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"

      Real*8 Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
      Call QEnter('SSB')
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
*---- SSB Exchange -- unlike OPTX, SSB has its LDA part included !
*
      Coeff=1.0d0

      Call xSSB(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          Coeff,iSpin,F_xc,T_X)
*
*---- PBE Correlation
*
      Coeff=1.0d0
      Call CPBE(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &         Coeff,iSpin,F_xc,T_X)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('SSB')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
