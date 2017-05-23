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
      Subroutine NDSD_Ts(mGrid,Rho,nRho,nDmat,F_xc,dF_dRho,
     &                   ndF_dRho,Coeff,T_X)
************************************************************************
*                                                                      *
* Object:  compute Func for Thomas-Fermi KE functional                 *
*          compute non-TF part (rho_B dependent) of NDSD potential     *
*                                                                      *
*          (see J.-M. Garcia Lastra, J. W. Kaminski, T. A. Wesolowski, *
*                                   J. Chem. Phys.  129 (2008) 074107.)*
*                                                                      *
*          Note: for a spin-polarized rho_B (environment density), the *
*                NDSD potential is computed using the alpha+beta       *
*                density, gradient and laplacian.                      *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
      Real*8 Fexp, Vt_lim
      External Fexp, Vt_lim
      Real*8 wGradRho(2:4)
*                                                                      *
************************************************************************
*                                                                      *
*
      Two3=Two/Three
      Five3=Five/Three
      Cf=(Three/Ten)*(three*Pi**Two)**Two3
      Rho_min=T_X*1.0D-2
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute value of energy and integrand on the grid
*                                                                      *
************************************************************************
*                                                                      *
      If (nDmat.eq.1) Then
         Do iGrid = 1, mGrid
            d_sys=Two*Rho(1,iGrid)
            If (d_sys.lt.T_X) Go To 100
*
*------- Kinetic energy contributions
*
            functional = Cf*d_sys**Five3
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*
*------- Contributions to the potential
*
            Do k=2,4
               wGradRho(k)=Two*Rho(k,iGrid)
            End Do
            wLaplRho=Two*Rho(6,iGrid)
*
            dfunc_NDSD = Fexp(d_sys,wGradRho(2))
     &                 * Vt_lim(d_sys,wGradRho(2),wLaplRho)
            dF_dRho(1,iGrid) = dF_dRho(1,iGrid)
     &                       + Coeff*dfunc_NDSD
*
 100        Continue
*
         End Do
*
      ElseIf (nDmat.eq.2) Then

         Cf = Cf*(Two**Two3)

         Do iGrid = 1, mGrid
            da_sys =Max(Rho_Min,Rho(1,iGrid))
            db_sys =Max(Rho_Min,Rho(2,iGrid))
            DTot=da_sys+db_sys
            If (DTot.lt.T_X) Go To 200
*
*------- Kinetic energy contributions
*
            functional=Cf*(da_sys**Five3+db_sys**Five3)
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*
*------- Contributions to the potential
*
            Do k=2,4
               ik=k+1
               jk=ik+3
               wGradRho(k)=Rho(ik,iGrid)+Rho(jk,iGrid)
            End Do
            wLaplRho=Rho(11,iGrid)+Rho(12,iGrid)
*
            dfunc_NDSD_alpha = Fexp(DTot,wGradRho(2))
     &                       * Vt_lim(DTot,wGradRho(2),wLaplRho)
            dfunc_NDSD_beta  = dfunc_NDSD_alpha
*
            dF_dRho(1,iGrid) = dF_dRho(1,iGrid)
     &                       + Coeff*dfunc_NDSD_alpha
            dF_dRho(2,iGrid) = dF_dRho(2,iGrid)
     &                       + Coeff*dfunc_NDSD_beta
*
 200        Continue
*
         End Do

      Else
         write(6,*) 'In NDSD_Ts: invalid # of densities. nDmat=  ',nDmat
         Call Abend()
      End If
*
      Return
      End
