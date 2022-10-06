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
!***********************************************************************
!                                                                      *
! This routine queries the existence of array data on runfile.         *
!                                                                      *
!***********************************************************************

subroutine LookUp_label(i,dtype,Label)

use RunFile_data, only: lw
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i
character(len=*) :: dtype
character(len=lw) :: Label
integer(kind=iwp) :: iTmp, nTmp
character(len=lw) :: RecLab(256)

!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call ffRun(dtype,nTmp,iTmp)
call cRdRun(dtype,RecLab,lw*256)
Label = RecLab(i)

return

end subroutine LookUp_label
