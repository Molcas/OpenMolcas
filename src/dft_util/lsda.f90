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
!***********************************************************************
#define _NEWCODE_
#ifdef _NEWCODE_
Subroutine LSDA(mGrid,nD)
use xc_f03_lib_m
use nq_Grid, only: F_xc => Exc, Rho
Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "hflda.fh"
#include "ksdft.fh"
Integer, parameter :: nFuncs=2
Real*8 :: Coeffs(nFuncs)=[1.0D0,1.0D0]
Integer*4, parameter :: func_id(nFuncs)=[int(1,4),int(8,4)]
TYPE(xc_f03_func_t) :: xc_func(nFuncs) ! xc functional
TYPE(xc_f03_func_info_t) :: xc_info(nFuncs) ! xc functional info
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Vosko-Wilk-Nusair correlation functional III
!
!---- Dirac exchange
!                                                                      *
!***********************************************************************
!                                                                      *
If (nD.eq.1) Rho(:,1:mGrid)=2.0D0*Rho(:,1:mGrid)
Do iFunc = 1, nFuncs
   ! Initialize libxc functional: nD = 2 means spin-polarized
   call xc_f03_func_init(xc_func(iFunc), func_id(iFunc), int(nD, 4))
   ! Get the functional's information
   xc_info(iFunc) = xc_f03_func_get_info(xc_func(iFunc))
End Do

Do iFunc = 1, nFuncs
   Coeff = Coeffs(iFunc)
   Select case(xc_f03_func_info_get_kind(xc_info(iFunc)))
     case (XC_EXCHANGE)
        Coeff = Coeff * CoefX
     case (XC_CORRELATION)
        Coeff = Coeff * CoefR
   End Select

   call libxc_interface(xc_func(iFunc),xc_info(iFunc),mGrid,nD,F_xc,Coeff)
End Do

Do iFunc = 1, nFuncs
   call xc_f03_func_end(xc_func(iFunc))
End Do
If (nD.eq.1) Rho(:,1:mGrid)=0.5D0*Rho(:,1:mGrid)
!                                                                      *
!***********************************************************************
!                                                                      *
Return
End Subroutine LSDA
#else
     Subroutine LSDA(mGrid,iSpin)
      use nq_Grid, only: F_xc => Exc
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Coeff=One*CoefR
      Call VWN_III(mGrid,iSpin,F_xc,Coeff)
      Coeff=One*CoefX
      Call DiracX(mGrid,iSpin,F_xc,Coeff)
      Return
      End
#endif

