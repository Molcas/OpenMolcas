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

subroutine Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Pren, Prem
integer(kind=iwp), intent(in) :: nBtch, mBtch, kBtch
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: nBtchV(3)
real(kind=wp) :: PrenV(2)

if (.not. Is_Real_Par()) return
if (nProcs == 1) return
PrenV(1) = Pren
PrenV(2) = Prem
call GADGOP(PrenV,2,'+')
Pren = PrenV(1)
Prem = PrenV(2)

nBtchV(1) = nBtch
nBtchV(2) = mBtch
nBtchV(3) = kBtch
call GAIGOP(nBtchV,3,'+')
nBtchV(1) = nBtch
nBtchV(2) = mBtch
nBtchV(3) = kBtch
#else
#include "macros.fh"
unused_var(Pren)
unused_var(Prem)
unused_var(nBtch)
unused_var(mBtch)
unused_var(kBtch)
#endif

end subroutine Sync_Data
