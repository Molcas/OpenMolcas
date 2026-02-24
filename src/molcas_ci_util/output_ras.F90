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

module OutPut_RAS
!------------------------------------------------------
! Used by any rasscf subroutine that writes output
!------------------------------------------------------
!***********************************************************************
!                                                                      *
!      Termination codes of the different program sections             *
!                                                                      *
!***********************************************************************
!
! Rc_CI   =  0 : CI-vectors are converged
!         = 16 : No convergence in the CI-section
!
! Rc_SX   =  0 : Super-CI method converged
!         = 16 : No convergence in the SX-section

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: IPRGLB, IPRLOC(7), Rc_CI, Rc_SX

public :: IPRGLB, IPRLOC, Rc_CI, Rc_SX

end module OutPut_RAS
