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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

function Fexp(rho,drho)
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!** Switching function of the NDSD potential:                        ***
!**     J.-M. Garcia Lastra, J. W. Kaminski, T. A. Wesolowski,       ***
!**                               J. Chem. Phys.  129 (2008) 074107. ***
!**                                                                  ***
!** Author: F. Aquilante, Geneva July 2010                           ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************

use Constants, only: One, Three, Two, Pi
use Definitions, only: wp

implicit none
real(kind=wp) :: Fexp
real(kind=wp), intent(in) :: rho, drho(3)
real(kind=wp) :: eir_rBmin, eis_sBmax, eis_sBmin, er_rBmin, es_sBmax, es_sBmin, fact, factinv, rhoinv, rhoinv13, sB, xnorm, xnorm_
real(kind=wp), parameter :: lambda = 5.0e2_wp, One3 = One/Three, rBmin = 0.7_wp, sBmin = 0.3_wp, sBmax = 0.9_wp

rhoinv = One/rho
rhoinv13 = rhoinv**One3
fact = Two*(Three*Pi**2)**One3
factinv = One/fact
xnorm = drho(1)**2+drho(2)**2+drho(3)**2
xnorm_ = sqrt(xnorm)
sB = factinv*rhoinv*rhoinv13*xnorm_

es_sBmin = exp(lambda*(sBmin-sB))
es_sBmax = exp(lambda*(sBmax-sB))
er_rBmin = exp(lambda*(rBmin-rho))

eis_sBmin = One/(es_sBmin+One)
eis_sBmax = One/(es_sBmax+One)
eir_rBmin = One/(er_rBmin+One)

Fexp = eis_sBmin*(One-eis_sBmax)*eir_rBmin

end function Fexp
