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
! Copyright (C) 2000-2016, Valera Veryazov                             *
!***********************************************************************
!*****************************************************************************
!                                                                            *
! Author:   Valera Veryazov 2000-2016                                        *
!           Theoretical Chemistry                                            *
!           Lund University                                                  *
!           Sweden                                                           *
!                                                                            *
!*****************************************************************************

subroutine systemf(c,rc)

implicit none
character*(*) C
character*1024 C2
integer LenC, StrnLn, i, rc

LenC = StrnLn(C)
if (LenC > 1024-1) then
  write(6,*) ' Error in systemf.f ! LenC :',LenC
  call abend()
end if

do i=1,lenc
  c2(i:i) = c(i:i)
end do
call systemc(c2,lenc,rc)

return

end subroutine systemf
