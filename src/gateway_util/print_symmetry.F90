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

use Symmetry_Info, only: iChTbl, iOper, lBsFnc, lIrrep, nIrrep, SymLab
use Definitions, only: iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: i, iIrrep, iPrint, iRout, j, jIrrep, nOper
character(len=80) :: frmt
character(len=*), parameter :: ChSymO(0:7) = ['  E  ', &
                                              's(yz)', &
                                              's(xz)', &
                                              'C2(z)', &
                                              's(xy)', &
                                              'C2(y)', &
                                              'C2(x)', &
                                              '  i  '], &
                               SymOpr(0:7) = [' Unit operation              ', &
                                              ' Reflection in the yz-plane  ', &
                                              ' Reflection in the xz-plane  ', &
                                              ' Rotation around the z-axis  ', &
                                              ' Reflection in the xy-plane  ', &
                                              ' Rotation around the y-axis  ', &
                                              ' Rotation around the x-axis  ', &
                                              ' Inversion through the origin']

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
!                                                                      *
!***********************************************************************
!                                                                      *
write(u6,*)
call CollapseOutput(1,'   Symmetry information:')
write(u6,'(3X,A)') '   ---------------------'
write(u6,*)
!                                                                      *
!***********************************************************************
!                                                                      *
if (nIrrep /= 1) then
  write(u6,'(19X,A)') ' --- Group Generators ---'
  nOper = 0
  if (nIrrep == 8) nOper = 3
  if (nIrrep == 4) nOper = 2
  if (nIrrep == 2) nOper = 1
  do i=1,nOper
    j = i
    if (i == 3) j = 4
    write(u6,'(19X,A)') SymOpr(iOper(j))
  end do
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
write(u6,'(19X,A,A)') ' Character Table for ',SymLab
write(u6,*)
write(frmt,'(A,I1,A)') '(20X,A3,1X,',nIrrep,'(1X,I5),2X,A)'
write(u6,'(27X,8(A5,1X))') (ChSymO(iOper(iIrrep)),iIrrep=0,nIrrep-1)
do iIrrep=0,nIrrep-1
  write(u6,frmt) lIrrep(iIrrep),(iChTbl(iIrrep,jIrrep),jIrrep=0,nIrrep-1),trim(lBsFnc(iIrrep))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call CollapseOutput(0,'  Symmetry information:')
write(u6,*)

return

end subroutine Print_Symmetry
