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
SUBROUTINE MkFock(CMO,nCMO,FIFA,nFIFA,DREF,nDREF,HONE,nHONE)
use caspt2_global, only: FIMO, FAMO
use caspt2_module, only: IfChol
use definitions, only: iwp, wp
IMPLICIT None
integer(kind=iwp), intent(in):: nCMO, nFIFA, nDREF, nHONE
real(kind=wp), intent(in):: CMO(nCMO), DREF(nDREF)
real(kind=wp), intent(inout):: FIFA(nFIFA)
real(kind=wp), intent(inout):: HONE(nHONE)

! Compute the Fock matrix in MO basis for state Jstate
! Fock matrix in MO basis: FIMO, FAMO, FIFA

if (IfChol) then
! INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis
   call INTCTL2(CMO,NCMO,DREF,nDREF,FIFA,nFIFA)
else
   CALL FMAT_CASPT2(FIMO,SIZE(FIMO),FAMO,SIZE(FAMO),DREF,nDREF,HONE)
end If

! Modify the Fock matrix if needed (G Family of modifications).
! You don't have to be beautiful to turn me on
CALL NEWFOCK(FIFA,nFIFA,CMO,NCMO,DREF,nDREF)

END SUBROUTINE MkFock
