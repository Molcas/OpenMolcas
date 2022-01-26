!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine NDSD_Ts(mGrid,nDmat)
!***********************************************************************
!                                                                      *
! Object:  compute Func for Thomas-Fermi KE functional                 *
!          compute non-TF part (rho_B dependent) of NDSD potential     *
!                                                                      *
!          (see J.-M. Garcia Lastra, J. W. Kaminski, T. A. Wesolowski, *
!                                   J. Chem. Phys.  129 (2008) 074107.)*
!                                                                      *
!          Note: for a spin-polarized rho_B (environment density), the *
!                NDSD potential is computed using the alpha+beta       *
!                density, gradient and laplacian.                      *
!                                                                      *
!***********************************************************************
      use nq_Grid, only: Rho, GradRho, Lapl
      use nq_Grid, only: vRho
      use nq_Grid, only: F_xc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Fexp, Vt_lim
      External Fexp, Vt_lim
      Real*8 wGradRho(1:3)
      Real*8, Parameter:: T_X=1.0D-20
      Real*8, Parameter:: Coeff=1.0D0
!                                                                      *
!***********************************************************************
!                                                                      *
      Two3=Two/Three
      Five3=Five/Three
      Cf=(Three/Ten)*(three*Pi**Two)**Two3
      Rho_min=T_X*1.0D-2
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Compute value of energy and integrand on the grid
!                                                                      *
!***********************************************************************
!                                                                      *
      If (nDmat.eq.1) Then
         Do iGrid = 1, mGrid
            d_sys=Two*Rho(1,iGrid)
            If (d_sys.lt.T_X) Go To 100
!
!------- Kinetic energy contributions
!
            functional = Cf*d_sys**Five3
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
!
!------- Contributions to the potential
!
            Do k=1,3
               wGradRho(k)=Two*GradRho(k,iGrid)
            End Do
            wLaplRho=Two*Lapl(1,iGrid)
!
            dfunc_NDSD = Fexp(d_sys,wGradRho(1))* Vt_lim(d_sys,wGradRho(1),wLaplRho)
            vRho(1,iGrid) = vRho(1,iGrid)+ Coeff*dfunc_NDSD
!
 100        Continue
!
         End Do
!
      ElseIf (nDmat.eq.2) Then

         Cf = Cf*(Two**Two3)

         Do iGrid = 1, mGrid
            da_sys =Max(Rho_Min,Rho(1,iGrid))
            db_sys =Max(Rho_Min,Rho(2,iGrid))
            DTot=da_sys+db_sys
            If (DTot.lt.T_X) Go To 200
!
!------- Kinetic energy contributions
!
            functional=Cf*(da_sys**Five3+db_sys**Five3)
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
!
!------- Contributions to the potential
!
            Do k=1,3
               wGradRho(k)=Rho(k,iGrid)+Rho(k+3,iGrid)
            End Do
            wLaplRho=Lapl(1,iGrid)+Lapl(2,iGrid)
!
            dfunc_NDSD_alpha = Fexp(DTot,wGradRho(1))* Vt_lim(DTot,wGradRho(1),wLaplRho)
            dfunc_NDSD_beta  = dfunc_NDSD_alpha
!
            vRho(1,iGrid) = vRho(1,iGrid)+ Coeff*dfunc_NDSD_alpha
            vRho(2,iGrid) = vRho(2,iGrid)+ Coeff*dfunc_NDSD_beta
!
 200        Continue
!
         End Do

      Else
         write(6,*) 'In NDSD_Ts: invalid # of densities. nDmat=  ',nDmat
         Call Abend()
      End If
!
      Return
      End
