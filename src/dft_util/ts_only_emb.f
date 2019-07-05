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
* Copyright (C) 2015, Francesco Aquilante                              *
*               2015, Alexander Zech                                   *
************************************************************************
      Subroutine Ts_only_emb(mGrid,Rho,nRho,P2_ontop,
     &                    nP2_ontop,nDmat,F_xc,
     &                    dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,
     &                    T_X)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "hflda.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       F_xc(mGrid),
     &       P2_ontop(nP2_ontop, mGrid),
     &       dF_dP2ontop(ndF_dP2ontop, mGrid)
      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
*
************************************************************************
*
*---- Thomas-Fermi Kinetic energy functional
*
      Coeff=One
      Call TF_Ts(mGrid,Rho,nRho,nDmat,F_xc,
     &                 dF_dRho,ndF_dRho,Coeff,T_X)
************************************************************************
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If

      Return
      End
