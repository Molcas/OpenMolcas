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
      Subroutine vW_Ts(mGrid,nDmat,F_xc,
     &                       Coeff,T_X)
************************************************************************
*                                                                      *
* Object:  compute Func for von Weizsacker KE functional               *
*          No potential computed!!!                                    *
*                                                                      *
************************************************************************
      use nq_grid, only: Rho, Sigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*
      One8=One/Eight
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
            snorm=Sigma(1,iGrid)
            functional = Half*snorm/d_sys
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*
*
 100        Continue
*
         End Do
*
      ElseIf (nDmat.eq.2) Then


         Do iGrid = 1, mGrid
            da_sys =Max(Rho_Min,Rho(1,iGrid))
            db_sys =Max(Rho_Min,Rho(2,iGrid))
            DTot=da_sys+db_sys
            If (DTot.lt.T_X) Go To 200
*
*------- Kinetic energy contributions
*
            snorm=Sigma(1,iGrid)
            functional = One8*snorm/da_sys
            snorm=Sigma(3,iGrid)
            functional = functional + One8*snorm/db_sys
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional

*
 200        Continue
*
         End Do

      Else
         write(6,*) 'In vW_Ts: invalid # of densities. nDmat=  ',nDmat
         Call Abend()
      End If
*
      Return
      End
