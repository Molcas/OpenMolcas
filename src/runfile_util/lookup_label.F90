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

subroutine LookUp_label(i,type,Label)

#include "pg_ca_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) type
character*16 Label
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
character*16 RecLab(256)
integer nTmp, iTmp
integer i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call ffRun(type,nTmp,iTmp)
call cRdRun(type,RecLab,16*256)
Label = RecLab(i)

return

end subroutine LookUp_label
