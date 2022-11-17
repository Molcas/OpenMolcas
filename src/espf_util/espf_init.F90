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

subroutine espf_init(natom,nAtQM,ipCord,ipIsMM,ipExt)
! ESPF initialization:
!   natom: number of atoms
!   ipCord: pointer to the atom coordinates
!   ipIsMM: pointer to the "logical" MM or QM array
!   ipExt: pointer to the external potential array

use espf_global, only: MxExtPotComp
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: natom, nAtQM, ipCord, ipIsMM, ipExt
#include "WrkSpc.fh"
integer(kind=iwp) :: nAtMM

call Get_iScalar('Unique atoms',natom)
call GetMem('AtomCoord','Allo','Real',ipCord,3*natom)
call Get_dArray('Unique Coordinates',Work(ipCord),3*natom)
call GetMem('IsMM for atoms','Allo','Inte',ipIsMM,natom)
call MMCount(natom,nAtMM,iWork(ipIsMM))
nAtQM = natom-nAtMM
call GetMem('ExtPot','ALLO','REAL',ipExt,natom*MxExtPotComp)
call dcopy_(MxExtPotComp*natom,[Zero],0,Work(ipExt),1)

return

end subroutine espf_init
