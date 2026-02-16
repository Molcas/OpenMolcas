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
!***********************************************************************
SUBROUTINE MkFock()
use caspt2_global, only: CMO, nCMO, FIFA
use caspt2_module, only: IfChol
use definitions, only: iwp, wp
IMPLICIT None

! Compute the Fock matrix in MO basis for state Jstate
! INTCTL1/INTCTL2 call TRACTL(0) and other routines to compute the
! Fock matrix in MO basis: FIMO, FAMO, FIFA and orbital energies

if (IfChol) then
! INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis
   call INTCTL2(CMO,NCMO)
else
! INTCTL1 uses TRAONE and FOCK_RPT2, to get the matrices in MO basis
   call INTCTL1(CMO,NCMO)
end If

! Modify the Fock matrix if needed (G Family of modifications).
! You don't have to be beautiful to turn me on
CALL NEWFOCK(FIFA,SIZE(FIFA),CMO,NCMO)

Call MkEPS()

END SUBROUTINE MkFock
