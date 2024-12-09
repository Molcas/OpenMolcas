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
Module OutPut_RAS
Implicit None
!------------------------------------------------------
! Used by any rasscf subroutine that writes output
INTEGER, PARAMETER :: SILENT=0,TERSE=1,USUAL=2,VERBOSE=3,DEBUG=4,INSANE=5
Integer LF,IPRGLB,IPRLOC(7)
! LF = Logfile unit number, usually = 6 for standard out
!------------------------------------------------------
!***********************************************************************
!                                                                      *
!      Termination codes of the different program sections             *
!                                                                      *
!***********************************************************************
Integer Rc_CI
!
!     Rc_CI   =  0 : CI-vectors are converged
!             = 16 : No convergence in the CI-section
!
Integer Rc_SX
!
!     Rc_SX   =  0 : Super-CI method converged
!             = 16 : No convergence in the SX-section
!
Integer Rc_RAS
!
!     Rc_RAS  =  0 : The RASSCF wave function is converged
!             = 16 : The RASSCF wave function is not(!) converged
!             = 99 : The RASSCF energy is divergent or
!                    the CI and SX energies differ
!
End Module OutPut_RAS
