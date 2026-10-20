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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine TRAFRO(MODE)

use caspt2_global, only: CMO, CMO_Internal, CMOPT2, NCMO
use general_data, only: nAsh
use caspt2_module, only: IfChol, NFRO, NISH, NORB, NOSH, NSSH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: MODE
integer(kind=iwp) :: nFroTmp(8), nOrbTmp(8), nOshTmp(8)

if (Mode == 1) then
  nFroTmp(1:nSym) = nFro(1:nSym)
  nOshTmp(1:nSym) = nOsh(1:nSym)
  nOrbTmp(1:nSym) = nOrb(1:nSym)
  nOsh(1:nSym) = nFro(1:nSym)+nIsh(1:nSym)+nAsh(1:nSym)
  nOrb(1:nSym) = nOsh(1:nSym)+nSsh(1:nSym)
  nFro(1:nSym) = 0
end if

call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
CMO => CMO_Internal
CMO(:) = CMOPT2(:)
if (IfChol) then
  call TRACHO3(CMO,NCMO)
else
  call TRACTL(nCMO,CMO,0)
end if
call mma_deallocate(CMO_Internal)
nullify(CMO)

if (Mode == 1) then
  nFro(1:nSym) = nFroTmp(1:nSym)
  nOsh(1:nSym) = nOshTmp(1:nSym)
  nOrb(1:nSym) = nOrbTmp(1:nSym)
end if

return

end subroutine TRAFRO
