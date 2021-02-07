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

subroutine ReadIn_Dynamix(Task,nTasks,mTasks)

implicit real*8(a-h,o-z)
integer Task(nTasks)

! Copy input from standard input to a local scratch file

LuSpool = isfreeunit(21)
call SpoolInp(LuSpool)

! Read input

#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix calls RdInp_Dynamix.'
#endif
call RdInp_Dynamix(LuSpool,Task,nTasks,mTasks)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix back from RdInp_Dynamix.'
#endif

! Remove local copy of standard input

close(LuSpool)

return

end subroutine ReadIn_Dynamix
