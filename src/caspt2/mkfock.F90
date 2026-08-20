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

subroutine MkFock(CMO,nCMO,FIMO,NFIMO,FIFA,nFIFA,DREF,nDREF,HONE,nHONE,INITIATE)

use caspt2_module, only: IfChol
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCMO, nFIMO, nFIFA, nDREF, nHONE
real(kind=wp), intent(in) :: CMO(nCMO), DREF(nDREF)
real(kind=wp), intent(inout) :: FIMO(nFIMO), FIFA(nFIFA), HONE(nHONE)
logical(kind=iwp), intent(inout) :: INITIATE

! Compute the Fock matrix in MO basis for state Jstate
! Fock matrix in MO basis: FIMO, FIFA

if (IfChol) then
  ! INTCTL2 uses TraCho2 to generate the fock matrix in AO basis. Subsequently, FMatCho
  ! transform to the MO basis.
  ! The one-electron Hamiltonian and G(D^{inactive+frozen}) are computed and added in TraCho2 (to FFIAO),
  ! so FIMO (and FIFA) remains unchanged.
  call INTCTL2(CMO,NCMO,DREF,nDREF,FIFA,nFIFA,FIMO,nFIMO)
  ! TraOne, which calculates HONE = one-electron Hamiltonian + G(D^frozen), is not necessary,
  ! because HONE is not used anywhere at present; if HONE is used in the future, this must be constructed.
  ! However, there should be little point in separating h + G(D^frozen) from G(D^inactive)
  if (Initiate) HONE(1:nHONE) = Zero
else
  ! Matrix elements generated directly from one-ham and two-electron integrals in th MO basis.
  if (Initiate) call TraOne(CMO,nCMO,HONE,nHONE)
  if (Initiate) call TraCtl(nCMO,CMO,0)
  call FMAT_CASPT2(FIFA,nFIFA,FIMO,nFIMO,DREF,nDREF,HONE,nHONE)
end if
INITIATE = .false.

! Modify the Fock matrix if needed (G Family of modifications).
! You don't have to be beautiful to turn me on
call NEWFOCK(FIFA,nFIFA,CMO,NCMO,DREF,nDREF)

end subroutine MkFock
