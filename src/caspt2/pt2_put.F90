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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine PT2_PUT(NSIZE,LAB,VEC)

use caspt2_global, only: LUDMAT
use caspt2_module, only: IADR10, cLab10
use definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp), intent(in) :: NSIZE
character(len=*), intent(in) :: LAB
real(kind=wp), intent(inout) :: VEC(NSIZE)
character(len=8) LAB1
integer(kind=iwp) I, IAD

I = 9-len(LAB)
if (I >= 1) then
  LAB1 = '        '
  LAB1(I:8) = LAB
else
  LAB1 = LAB(1:8)
end if

! FIND DISK ADDRESS:
do I=1,size(CLAB10)
  if (CLAB10(I) == '   EMPTY') then
    CLAB10(I) = LAB1
    IAD = IADR10(I,1)
    IADR10(I,2) = NSIZE
    call DDAFILE(LUDMAT,1,VEC,NSIZE,IAD)
    if (I < size(CLAB10)) IADR10(I+1,1) = IAD
    return
  else if (CLAB10(I) == LAB1) then
    if (NSIZE > IADR10(I,2)) then
      write(u6,*) ' ATTEMPT TO INCREASE SIZE OF A FIELD.'
      write(u6,*) ' SUBROUTINE PUT FAILS.'
      call ABEND()
    end if
    IAD = IADR10(I,1)
    IADR10(I,2) = NSIZE
    call DDAFILE(LUDMAT,1,VEC,NSIZE,IAD)
    return
  end if
end do

write(u6,*) ' NO MORE AVAILABLE FIELDS ON FILE DENS.'
write(u6,*) ' SUBROUTINE PUT FAILS.'
call ABEND()

end subroutine PT2_PUT
