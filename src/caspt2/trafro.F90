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
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NASH, NOSH, NSSH, NORB
use definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: MODE
integer(kind=iwp) :: nFroTmp(8), nOshTmp(8), nOrbTmp(8)
integer(kind=iwp) :: jSym

if (Mode == 1) then
  do jSym=1,nSym
    nFroTmp(jSym) = nFro(jSym)
    nOshTmp(jSym) = nOsh(jSym)
    nOrbTmp(jSym) = nOrb(jSym)
    nOsh(jSym) = nFro(jSym)+nIsh(jSym)+nAsh(jSym)
    nOrb(jSym) = nOsh(jSym)+nSsh(jSym)
    nFro(jSym) = 0
  end do
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
  do jSym=1,nSym
    nFro(jSym) = nFroTmp(jSym)
    nOsh(jSym) = nOshTmp(jSym)
    nOrb(jSym) = nOrbTmp(jSym)
  end do
end if

return

end subroutine TRAFRO
