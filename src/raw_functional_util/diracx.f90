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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
      Subroutine DiracX(mGrid,nD,F_xc,Coeff)
      use xc_f03_lib_m
      use nq_Grid, only: Rho, vRho, l_casdft
      use KSDFT_Info, only: F_xca, F_xcb
      implicit none
      Real*8 :: F_xc(mGrid)
      Real*8 :: Coeff
      integer :: mgrid, nD, iGrid

      ! Work memory for libxc
      Real*8 :: func(mGrid), dfunc_drho(nD,mGrid)

      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! Slater exchange
      integer*4, parameter :: func_id = 1

      ! Initialize memory
      func(:) = 0.0
      dfunc_drho(:,:) = 0.0

      ! Initialize libxc functional: nRho = 2 means spin-polarized
      call xc_f03_func_init(xc_func, func_id, int(nRho, 4))

      ! Get the functional's information
      xc_info = xc_f03_func_get_info(xc_func)

      ! Evaluate energy depending on the family
      select case (xc_f03_func_info_get_family(xc_info))
      case(XC_FAMILY_LDA)
         call xc_f03_lda_exc_vxc(xc_func, mGrid, Rho(1,1), func(1), dfunc_drho(1,1))
      end select

      ! Libxc evaluates energy density per particle; multiply by
      ! density to get out what we really want
      ! Collect the potential
      If (nD.eq.1) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)*Rho(1, iGrid)
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)*(Rho(1, iGrid) + Rho(2, iGrid))
            vRho(1,iGrid) = vRho(1,iGrid) + Coeff*dfunc_drho(1, iGrid)
            vRho(2,iGrid) = vRho(2,iGrid) + Coeff*dfunc_drho(2, iGrid)
         End Do
         If (l_casdft) Then
            ! more to come here.
         End If
      End If

      call xc_f03_func_end(xc_func)
      Return

    End Subroutine DiracX
