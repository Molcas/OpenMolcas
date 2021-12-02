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
      Subroutine DiracX(mGrid,Rho,nRho,iSpin,F_xc,dF_dRho,
     &                  ndF_dRho,Coeff,T_X)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      Use xc_f03_lib_m
      Use Definitions, Only: LibxcInt, LibxcReal, LibxcSize
      Implicit None
      Real*8 :: Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
      Real*8 :: Coeff, T_X
      Integer :: mGrid, nrho, ndf_drho, ispin, iGrid
*     Work memory for libxc
      Real(kind=LibxcReal) :: func(mGrid), dfunc_drho(iSpin,mGrid)
*     xc functional
      Type(xc_f03_func_t) :: xc_func
*     xc functional info
      Type(xc_f03_func_info_t) :: xc_info
*     Slater exchange
      Integer(kind=LibxcInt), Parameter :: func_id=1
      Integer, Parameter :: ipR=1, ipRa=1, ipRb=2
#include "real.fh"
*
*---- Initialize memory
*
      func(:) = Zero
      dfunc_drho(:,:) = Zero
*
*---- Initialize libxc functional: nRho = 2 means spin-polarized
*
      Call xc_f03_func_init(xc_func, func_id, Int(nRho, LibxcInt))
*
*---- Get the functional's information
*
      xc_info = xc_f03_func_get_info(xc_func)
*
*---- Evaluate energy depending on the family
*
      Select Case (xc_f03_func_info_get_family(xc_info))
         Case(XC_FAMILY_LDA)
            Call xc_f03_lda_exc_vxc(xc_func, Int(mGrid,kind=LibxcSize),
     &                              Real(Rho(:,:),kind=LibxcReal),
     &                              func(:), dfunc_drho(:,:))
      End Select
*
*---- Libxc evaluates energy density per particle;
*     multiply by density to get out what we really want
*
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            func(iGrid) = func(iGrid) * Rho(1, iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            func(iGrid) = func(iGrid) * (Rho(1, iGrid) + Rho(2, iGrid))
         End Do
      End If
*
*---- Collect the potential
*
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)
            dF_dRho(ipR,iGrid) = dF_dRho(ipR,iGrid) +
     &                           Coeff*dfunc_drho(1, iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)
            dF_dRho(ipRa,iGrid) = dF_dRho(ipRa,iGrid) +
     &                            Coeff*dfunc_drho(1, iGrid)
            dF_dRho(ipRb,iGrid) = dF_dRho(ipRb,iGrid) +
     &                            Coeff*dfunc_drho(2, iGrid)
         End Do
      End If
*
      Return
* Avoid unused argument warnings
      If (.False.) Call Unused_real(T_X)
      End
