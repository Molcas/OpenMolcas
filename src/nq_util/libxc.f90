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
Module libxc
Real*8, Allocatable:: func(:)
Real*8, Allocatable:: dfunc_drho(:,:)
Real*8, Allocatable:: dfunc_dsigma(:,:)
Real*8, Allocatable:: dfunc_dtau(:,:)
Real*8, Allocatable:: dfunc_dlapl(:,:)
Logical :: Only_exc=.False.
End Module libxc
