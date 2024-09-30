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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************

subroutine Mk_EOrb()

use InfSCF, only: CMO, EOrb, FockAO, nBas, nD, nOrb, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iD, nCMO, nEOrb, nFck

nFck = size(FockAO,1)
nEOrb = size(EOrb,1)
nCMO = size(CMO,1)

do iD=1,nD
  call MkEorb_Inner(FockAO(:,iD),nFck,CMO(:,iD),nCMO,EOrb(:,iD),nEorb,nSym,nBas,nOrb)
  if (iD == 1) then
    call Put_darray('OrbE',Eorb(:,iD),size(EOrb,1))
  else
    call Put_darray('OrbE_ab',Eorb(:,iD),size(EOrb,1))
  end if
end do

return

end subroutine Mk_EOrb
