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
! Copyright (C) Jeppe Olsen                                            *
!***********************************************************************

subroutine BLKCLS(IBLKS,NBLKS,IBLKCLS,ISPSPCL,NCLS,LCLS,NOCTPA,NOCTPB)
! Class of each block, and dimension of each class
!
! Jeppe Olsen

implicit real*8(A-H,O-Z)
! Input
integer IBLKS(8,NBLKS)
integer ISPSPCL(NOCTPA,NOCTPB)
! Output
integer IBLKCLS(NBLKS), LCLS(NCLS)

!write(6,*) ' ISPSPCL'
!call IWRTMA(ISPSPCL,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
IZERO = 0
call ISETVC(LCLS,IZERO,NCLS)
do JBLK=1,NBLKS
  IICLS = ISPSPCL(IBLKS(1,JBLK),IBLKS(2,JBLK))
  IBLKCLS(JBLK) = IICLS
  LCLS(IICLS) = LCLS(IICLS)+IBLKS(8,JBLK)
end do

NTEST = 0
if (NTEST >= 100) then
  write(6,*)
  write(6,*) ' BLKCLS Speaking'
  write(6,*) ' ==============='
  write(6,*)
  write(6,*) ' Dimension of each class'
  call IWRTMA(LCLS,1,NCLS,1,NCLS)
  write(6,*)
  write(6,*) ' Class of each block :'
  call IWRTMA(IBLKCLS,1,NBLKS,1,NBLKS)
end if

end subroutine BLKCLS
