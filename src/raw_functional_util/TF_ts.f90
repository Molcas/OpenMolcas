!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
#define _NEWCODE_
#ifdef _NEWCODE_
      Subroutine TF_Ts(mGrid,nD,F_xc,Coeff)
      use xc_f03_lib_m
      use nq_Grid, only: Rho
      use libxc
      implicit none
      integer :: mGrid, nD, nRho
      Real*8 :: F_xc(mGrid)
      Real*8 :: Coeff

      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! Slater exchange
      integer*4, parameter :: func_id = 50

      nRho=SIZE(Rho,1)
      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      If (nD.eq.1) Rho(:,1:mGrid)=2.0D0*Rho(:,1:mGrid)

      call libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

      call xc_f03_func_end(xc_func)

      If (nD.eq.1) Rho(:,1:mGrid)=0.5D0*Rho(:,1:mGrid)
!     Call RecPrt('F_xc',' ',F_xc,1,mGrid)
!     Stop 123
      Return

    End Subroutine TF_Ts
#else
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
      Subroutine TF_Ts(mGrid,nDmat,F_xc,Coeff)
!***********************************************************************
!                                                                      *
! Object:  compute Func and potential for Thomas-Fermi KE functional   *
!                                                                      *
!***********************************************************************
      use nq_Grid, only: Rho
      use nq_Grid, only: vRho
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20
!                                                                      *
!***********************************************************************
!                                                                      *
!
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
!           functional = Cf*d_sys**Five3
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
!
!------- Contributions to the potential
!
            dfunc_TF = Five3*Cf*d_sys**Two3
            vRho(1,iGrid) = vRho(1,iGrid)+ Coeff*dfunc_TF
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
            dfunc_TF_alpha=Five3*Cf*da_sys**Two3
            dfunc_TF_beta =Five3*Cf*db_sys**Two3
            vRho(1,iGrid) = vRho(1,iGrid)+ Coeff*dfunc_TF_alpha
            vRho(2,iGrid) = vRho(2,iGrid)+ Coeff*dfunc_TF_beta
!
 200        Continue
!
         End Do

      Else
         write(6,*) 'In TF_Ts: invalid # of densities. nDmat=  ',nDmat
         Call Abend()
      End If
!
      Return
      End
#endif
