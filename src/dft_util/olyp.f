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
* Copyright (C) 2009, Grigory A. Shamov                                *
************************************************************************
      Subroutine OLYP(mGrid,Rho,nRho,
     &                iSpin,F_xc)
************************************************************************
*                                                                      *
* Object:    OLYP combination                                          *
*            ref: Handy, Cohen, J. Mol.Phys 99 (2001) 403              *
*                                                                      *
* Author:    Grigory A Shamov, University of Manitoba 2009             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),
     &       F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*---- Dirac Exchange with OLYP factor!
*
      Coeff= 1.051510d0*CoefX
      Call Diracx(mGrid,iSpin,F_xc,Coeff)
*
*---- OPTX Exchange
*
      Coeff= 1.431690d0*CoefX
      Call xOPT(mGrid,
     &          Coeff,iSpin,F_xc)
*
*---- Lee-Yang-Parr Correlation
*
      Coeff=One*CoefR
      Call LYP(mGrid,
     &         Coeff,iSpin,F_xc)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_Integer(nRho)
         Call Unused_real_array(Rho)
      End If
      End
