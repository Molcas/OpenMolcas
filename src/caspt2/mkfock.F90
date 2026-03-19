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
SUBROUTINE MkFock(CMO,nCMO,FIMO,NFIMO,FIFA,nFIFA,DREF,nDREF,HONE,nHONE,INITIATE)
use caspt2_module, only: IfChol
use definitions, only: iwp, wp
IMPLICIT None
integer(kind=iwp), intent(in):: nCMO, nFIMO, nFIFA, nDREF, nHONE
real(kind=wp), intent(in):: CMO(nCMO), DREF(nDREF)
real(kind=wp), intent(inout):: FIMO(nFIMO), FIFA(nFIFA)
real(kind=wp), intent(inout):: HONE(nHONE)
logical(kind=iwp), intent(inout):: INITIATE

If (INITIATE) Call TraOne(CMO,nCMO,HONE,nHONE)

! Compute the Fock matrix in MO basis for state Jstate
! Fock matrix in MO basis: FIMO, FIFA

if (IfChol) then
!  INTCTL2 uses TraCho2 to generate the fock matrix in AO basis. Subsequently, FMatCho
!  transform to the MO basis.
   call INTCTL2(CMO,NCMO,DREF,nDREF,FIFA,nFIFA,HONE,nHONE,FIMO,nFIMO)
else
!  Matrix elements generated directly from one-ham and two-electron integrals in th MO basis.
   If (Initiate) Call TraCtl(nCMO,CMO,0)
   CALL FMAT_CASPT2(FIFA,nFIFA,FIMO,nFIMO,DREF,nDREF,HONE,nHONE)
end If
INITIATE=.FALSE.

! Modify the Fock matrix if needed (G Family of modifications).
! You don't have to be beautiful to turn me on
CALL NEWFOCK(FIFA,nFIFA,CMO,NCMO,DREF,nDREF)

END SUBROUTINE MkFock
