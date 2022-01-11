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
Subroutine LSDA(mGrid,nD)
use nq_Grid, only: F_xc => Exc, Rho
use libxc_parameters
Implicit None
Integer mGrid,nD
Real*8 Coeff
#include "real.fh"

! Specify DFT functionals and their corresponding coefficients
nFuncs=2
func_id(1:nFuncs)=[int(1,4),int(8,4)]
Coeffs(1:nFuncs)=[1.0D0,1.0D0]

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

Call Initiate_Libxc_functionals(nD)

Do iFunc = 1, nFuncs
   Coeff = Coeffs(iFunc)
   call libxc_interface(xc_func(iFunc),xc_info(iFunc),mGrid,nD,F_xc,Coeff)
End Do

Call Remove_Libxc_functionals()

If (nD.eq.1) Rho(:,1:mGrid)=0.5D0*Rho(:,1:mGrid)
!                                                                      *
!***********************************************************************
!                                                                      *
Return
End Subroutine LSDA
