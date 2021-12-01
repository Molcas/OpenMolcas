! ************************************************************************
! * This file is part of OpenMolcas.                                     *
! *                                                                      *
! * OpenMolcas is free software; you can redistribute it and/or modify   *
! * it under the terms of the GNU Lesser General Public License, v. 2.1. *
! * OpenMolcas is distributed in the hope that it will be useful, but it *
! * is provided "as is" and without any express or implied warranties.   *
! * For more details see the full text of the license in the file        *
! * LICENSE or in <http://www.gnu.org/licenses/>.                        *
! *                                                                      *
! * Copyright (C) 2000, Roland Lindh                                     *
! ************************************************************************
      Subroutine DiracX(mGrid,Rho,nRho,iSpin,F_xc,dF_dRho,ndF_dRho,Coeff,T_X)
! ************************************************************************
! *      Author:Roland Lindh, Department of Chemical Physics, University *
! *             of Lund, SWEDEN. November 2000                           *
! ************************************************************************
      use KSDFT_Info, only: F_xca, F_xcb
      use xc_f03_lib_m
      implicit none
!#include "real.fh"
!#include "nq_index.fh"
!#include "ksdft.fh"
      Real*8 :: Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
      Real*8 :: Coeff, T_X
      integer :: mgrid, nrho, ndf_drho, ispin, iGrid
      ! Work memory for libxc
      Real*8 :: func(mGrid), dfunc_drho(iSpin,mGrid)
      ! xc functional
      TYPE(xc_f03_func_t) :: xc_func
      ! xc functional info
      TYPE(xc_f03_func_info_t) :: xc_info

      ! Slater exchange
      integer*4, parameter :: func_id = 1
      integer, parameter :: ipR=1, ipRa=1, ipRb=2

      ! Initialize memory
      func = 0.0
      dfunc_drho = 0.0

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
      if (iSpin.eq.1) then
         do iGrid = 1, mGrid
            func(iGrid) = func(iGrid) * Rho(1, iGrid)
         end do
      else
         do iGrid = 1, mGrid
            func(iGrid) = func(iGrid) * (Rho(1, iGrid) + Rho(2, iGrid))
         end do
      end if

      ! Collect the potential
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            F_xc(iGrid) = F_xc(iGrid) + Coeff*func(iGrid)
            dF_dRho(ipR,iGrid) = dF_dRho(ipR,iGrid) + Coeff*dfunc_drho(1, iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            F_xc(iGrid) =F_xc(iGrid) +Coeff*func(iGrid)
            dF_dRho(ipRa,iGrid) = dF_dRho(ipRa,iGrid) + Coeff*dfunc_drho(1, iGrid)
            dF_dRho(ipRb,iGrid) = dF_dRho(ipRb,iGrid) + Coeff*dfunc_drho(2, iGrid)
         End Do
      End If
      Return
    End Subroutine DiracX
