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
! Copyright (C) 2000,2022, Roland Lindh                                *
!               2022, Susi Lehtola                                     *
!***********************************************************************
Module libxc_parameters
use xc_f03_lib_m
use Definitions, only: LibxcInt
Implicit None
#include "ksdft.fh"

Integer, parameter :: nFuncs_max=4
Integer :: i
Integer :: nFuncs=0
Real*8 :: Coeffs(nFuncs_Max)=[(0.0D0,i=1,nFuncs_Max)]
Integer(kind=LibxcInt) :: func_id(nFuncs_Max)=[(0_LibxcInt,i=1,nFuncs_Max)]

TYPE(xc_f03_func_t)      :: xc_func(nFuncs_Max) ! xc functional
TYPE(xc_f03_func_info_t) :: xc_info(nFuncs_Max) ! xc functional info

!
!***********************************************************************
!
Contains
!
!***********************************************************************
!
Subroutine Initiate_Libxc_functionals(nD)
use nq_Grid, only: l_casdft
Implicit None
Integer nD, iFunc
Real*8 :: Coeff

! if it is a mixed functional and we do MC-PDFT split it up in the components for
! further analysis.
If (nFuncs==1 .and. l_casdft) Then
   call xc_f03_func_init(xc_func(1), func_id(1), int(nD, kind=LibxcInt))
   nFuncs = Max(1,INT(xc_f03_num_aux_funcs(xc_func(1))))

   If (nFuncs/=1) Then
      call xc_f03_aux_func_ids(xc_func(1), func_id)
      call xc_f03_aux_func_weights(xc_func(1), Coeffs)
      call xc_f03_func_end(xc_func(1))
   End If

End If
Do iFunc = 1, nFuncs
   ! Initialize libxc functional: nD = 2 means spin-polarized
   call xc_f03_func_init(xc_func(iFunc), func_id(iFunc), int(nD, kind=LibxcInt))
   ! Get the functional's information
   xc_info(iFunc) = xc_f03_func_get_info(xc_func(iFunc))

! Reset coefficiants according to input

   Coeff = Coeffs(iFunc)
   Select case(xc_f03_func_info_get_kind(xc_info(iFunc)))
     case (XC_EXCHANGE)
        Coeff = Coeff * CoefX
     case (XC_CORRELATION)
        Coeff = Coeff * CoefR
   End Select
   Coeffs(iFunc) = Coeff

End Do

End Subroutine Initiate_Libxc_functionals
!
!***********************************************************************
!
Subroutine Remove_Libxc_functionals()
Implicit None
Integer iFunc
Do iFunc = 1, nFuncs
   call xc_f03_func_end(xc_func(iFunc))
End Do
Coeffs(:)=0.0D0
func_id(:)=0
End Subroutine Remove_Libxc_functionals
!
!***********************************************************************
!
Subroutine libxc_functionals(mGrid,nD)
use nq_Grid, only: F_xc
use nq_Grid, only: Rho, Sigma, Tau, Lapl
use nq_Grid, only:     vSigma,vTau
Implicit None
Integer mGrid,nD, iFunc
Real*8 Coeff

If (nD.eq.1) Then
                          Rho(:,1:mGrid)   =2.00D0*Rho(:,1:mGrid)
   If (Allocated(Sigma))  Sigma(:,1:mGrid) =4.00D0*Sigma(:,1:mGrid)
   If (Allocated(vSigma)) vSigma(:,1:mGrid)=0.50D0*vSigma(:,1:mGrid)
   If (Allocated(Lapl))   Lapl(:,1:mGrid)  =2.00D0*Lapl(:,1:mGrid)
   If (Allocated(vTau))   vTau(:,1:mGrid)  =2.00D0*vTau(:,1:mGrid)
Else
   If (Allocated(Tau))    Tau(:,1:mGrid)   =0.50D0*Tau(:,1:mGrid)
   If (Allocated(vTau))   vTau(:,1:mGrid)  =2.00D0*vTau(:,1:mGrid)
End If

Do iFunc = 1, nFuncs
   Coeff = Coeffs(iFunc)
   call libxc_interface(xc_func(iFunc),xc_info(iFunc),mGrid,nD,F_xc,Coeff)
End Do

If (nD.eq.1) Then
                          Rho(:,1:mGrid)   =0.50D0*Rho(:,1:mGrid)
   If (Allocated(Sigma))  Sigma(:,1:mGrid) =0.25D0*Sigma(:,1:mGrid)
   If (Allocated(vSigma)) vSigma(:,1:mGrid)=2.00D0*vSigma(:,1:mGrid)
   If (Allocated(Lapl))   Lapl(:,1:mGrid)  =0.50D0*Lapl(:,1:mGrid)
   If (Allocated(vTau))   vTau(:,1:mGrid)  =0.50D0*vTau(:,1:mGrid)
Else
   If (Allocated(Tau))    Tau(:,1:mGrid)   =2.00D0*Tau(:,1:mGrid)
   If (Allocated(vTau))   vTau(:,1:mGrid)  =0.50D0*vTau(:,1:mGrid)
End If

Return
End Subroutine libxc_functionals
!
!***********************************************************************
!
End Module libxc_parameters
