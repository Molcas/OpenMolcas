#define _NEWCODE_
#ifdef _NEWCODE_
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
     Subroutine vW_Ts(mGrid,Coeff,nD,F_xc)
      use xc_f03_lib_m
      use nq_Grid, only: Rho, Sigma
      use nq_Grid, only: vSigma
      use libxc
      implicit none
      integer :: mGrid, nD, nRho
      Real*8 :: F_xc(mGrid)
      Real*8 :: Coeff

      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! Becke exchange 88 exchange
      integer*4, parameter :: func_id = 52

      nRho=SIZE(Rho,1)
      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=2.00D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=4.00D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=0.50D0*vSigma(:,1:mGrid)
      End If

      call libxc_interface(xc_func,xc_info,mGrid,nD,F_xc,Coeff)

      call xc_f03_func_end(xc_func)

      If (nD.eq.1) Then
         Rho(:,1:mGrid)=0.50D0*Rho(:,1:mGrid)
         Sigma(:,1:mGrid)=0.25D0*Sigma(:,1:mGrid)
         vSigma(:,1:mGrid)=2.00D0*vSigma(:,1:mGrid)
      End If

    End Subroutine vW_Ts
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
      Subroutine vW_Ts(mGrid,nDmat,F_xc,Coeff)
!***********************************************************************
!                                                                      *
! Object:  compute Func for von Weizsacker KE functional               *
!          No potential computed!!!                                    *
!                                                                      *
!***********************************************************************
      use nq_grid, only: Rho, Sigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20
!                                                                      *
!***********************************************************************
!                                                                      *
!
      One8=One/Eight
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
            snorm=Sigma(1,iGrid)
            functional = Half*snorm/d_sys
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
!
!
 100        Continue
!
         End Do
!
      ElseIf (nDmat.eq.2) Then


         Do iGrid = 1, mGrid
            da_sys =Max(Rho_Min,Rho(1,iGrid))
            db_sys =Max(Rho_Min,Rho(2,iGrid))
            DTot=da_sys+db_sys
            If (DTot.lt.T_X) Go To 200
!
!------- Kinetic energy contributions
!
            snorm=Sigma(1,iGrid)
            functional = One8*snorm/da_sys
            snorm=Sigma(3,iGrid)
            functional = functional + One8*snorm/db_sys
            F_xc(iGrid)=F_xc(iGrid)+Coeff*functional

!
 200        Continue
!
         End Do

      Else
         write(6,*) 'In vW_Ts: invalid # of densities. nDmat=  ',nDmat
         Call Abend()
      End If
!
      Return
      End
#endif
