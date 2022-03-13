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

subroutine Memory_PolPrep(Que,ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri,nPol,nPart)

use Constants, only: Zero
use Definitions, only: iwp

implicit none
character(len=*) :: Que
integer(kind=iwp) :: ixx, iyy, izz, irr3, ixxi, iyyi, izzi, iGri, nPol, nPart
#include "WrkSpc.fh"
integer(kind=iwp) :: nSize

nSize = nPart*nPol
call GetMem('xx',Que,'Real',ixx,nSize**2)
call GetMem('yy',Que,'Real',iyy,nSize**2)
call GetMem('zz',Que,'Real',izz,nSize**2)
call GetMem('ixx',Que,'Real',ixxi,nSize**2)
call GetMem('iyy',Que,'Real',iyyi,nSize**2)
call GetMem('izz',Que,'Real',izzi,nSize**2)
call GetMem('irr3',Que,'Real',irr3,nSize**2)
call GetMem('iGri',Que,'Real',iGri,nSize**2)

if (Que(1:4) == 'Allo') then
  call dcopy_(nSize**2,[Zero],0,Work(ixx),1)
  call dcopy_(nSize**2,[Zero],0,Work(iyy),1)
  call dcopy_(nSize**2,[Zero],0,Work(izz),1)
  call dcopy_(nSize**2,[Zero],0,Work(ixxi),1)
  call dcopy_(nSize**2,[Zero],0,Work(iyyi),1)
  call dcopy_(nSize**2,[Zero],0,Work(izzi),1)
  call dcopy_(nSize**2,[Zero],0,Work(irr3),1)
  call dcopy_(nSize**2,[Zero],0,Work(iGri),1)
end if

return

end subroutine Memory_PolPrep
