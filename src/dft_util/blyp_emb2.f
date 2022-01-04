***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine BLYP_emb2(mGrid,
     &                     nDmat,F_xc)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
************************************************************************
      use OFembed, only: KEonly
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "hflda.fh"
      Real*8 F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*---- Thomas-Fermi Kinetic energy functional
*---- NDSD potential
*
      Coeff=One
      Call ndsd_Ts(mGrid,nDmat,F_xc,Coeff)

      If (KEonly) Return
*
*---- BLYP for exchange-correlation energy functional
*
      Call BLYP_(mGrid,nDmat,F_xc)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine BLYP_(mGrid,iSpin,F_xc)

      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange
*
      Coeff=Zero
      Call Diracx_ofe(mGrid,iSpin,F_xc,Coeff)
*                                                                      *
*---- Becke 88 Exchange
*
      Coeff=Zero
      Call xB88_ofe(mGrid,
     &              Coeff,iSpin,F_xc)
*
*---- Lee-Yang-Parr Correlation
*
      Coeff=Zero
      Call LYP_ofe(mGrid,
     &             Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
