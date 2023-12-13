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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine PXMem( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Integral_interfaces, only: int_mem
use Property_Label, only: PLabel
use Definitions, only: iwp, u6

implicit none
#include "mem_interface.fh"
procedure(int_mem) :: CntMem, EFMem, MltMem, NAMem

if (PLabel == 'NAInt ') then
  call PVMem(nHer,Mem,la,lb,lr,NAMem)
else if (PLabel == 'MltInt') then
  call PVMem(nHer,Mem,la,lb,lr,MltMem)
else if (PLabel == 'EFInt ') then
  call PVMem(nHer,Mem,la,lb,lr,EFMem)
else if (PLabel == 'CntInt') then
  call PVMem(nHer,Mem,la,lb,lr,CntMem)
else
  call WarningMessage(2,'PXMem: Illegal type!')
  write(u6,*) '       PLabel=',PLabel
  call Abend()
end if

return

end subroutine PXMem
