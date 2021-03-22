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

subroutine RPA_SetInc()

use Constants, only: Zero
use Definitions, only: iwp

implicit none
#include "rpa_config.fh"
#include "rpa_data.fh"
integer(kind=iwp) :: i, j

! rpa_config
Reference = 'Non'
RPAModel = 'None@Non'
DFTFunctional = 'Not defined     '
dRPA = .false.
SOSEX = .false.
doCD = .false.
doDF = .false.
doLDF = .false.
LumOrb = .false.
iPrint = 0
! rpa_data
do i=1,mTitle
  write(Title(i),'(80A1)')(' ',j=1,80)
end do
nTitle = 0
nSym = 0
call iZero(nFreeze,2)
call iZero(nBas,8)
call iZero(nOrb,8)
call iZero(nFro,16)
call iZero(nDel,16)
call iZero(nOcc,16)
call iZero(nVir,16)
call iZero(ip_CMO,2)
call iZero(l_CMO,2)
call iZero(ip_EMO,2)
call iZero(l_EMO,2)
call iZero(ip_OccEn,2)
call iZero(l_OccEn,2)
call iZero(ip_VirEn,2)
call iZero(l_VirEn,2)
NuclearRepulsionEnergy(1) = Zero

end subroutine RPA_SetInc
