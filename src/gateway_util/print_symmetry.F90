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
! Copyright (C) 2006, Roland Lindh                                     *
!***********************************************************************

subroutine Print_Symmetry()
!***********************************************************************
!                                                                      *
!     Object: to write the output of seward                            *
!                                                                      *
!     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
!             September '06                                            *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl, iOper, lIrrep, lBsFnc,SymLab

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
character SymOpr(0:7)*29, format*80, ChSymO(0:7)*5
data SymOpr/' Unit operation              ',' Reflection in the yz-plane  ',' Reflection in the xz-plane  ', &
            ' Rotation around the z-axis  ',' Reflection in the xy-plane  ',' Rotation around the y-axis  ', &
            ' Rotation around the x-axis  ',' Inversion through the origin'/
data ChSymO/'  E  ','s(yz)','s(xz)','C2(z)','s(xy)','C2(y)','C2(x)','  i  '/

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
LuWr = 6
!                                                                      *
!***********************************************************************
!                                                                      *
write(LuWr,*)
call CollapseOutput(1,'   Symmetry information:')
write(LuWr,'(3X,A)') '   ---------------------'
write(LuWr,*)
!                                                                      *
!***********************************************************************
!                                                                      *
if (nIrrep /= 1) then
  write(LuWr,'(19X,A)') ' --- Group Generators ---'
  nOper = 0
  if (nIrrep == 8) nOper = 3
  if (nIrrep == 4) nOper = 2
  if (nIrrep == 2) nOper = 1
  do i=1,nOper
    j = i
    if (i == 3) j = 4
    write(LuWr,'(19X,A)') SymOpr(iOper(j))
  end do
  write(LuWr,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
write(LuWr,'(19X,A,A)') ' Character Table for ',SymLab
write(LuWr,*)
write(format,'(A,I1,A)') '(20X,A3,1X,',nIrrep,'(1X,I5),2X,A)'
write(LuWr,'(27X,8(A5,1X))') (ChSymO(iOper(iIrrep)),iIrrep=0,nIrrep-1)
do iIrrep=0,nIrrep-1
  LenlBs = len(lBsFnc(iIrrep))
  write(LuWr,format) lIrrep(iIrrep),(iChTbl(iIrrep,jIrrep),jIrrep=0,nIrrep-1),lBsFnc(iIrrep)(1:iCLast(lBsFnc(iIrrep),LenlBs))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call CollapseOutput(0,'  Symmetry information:')
write(LuWr,*)

return

end subroutine Print_Symmetry
