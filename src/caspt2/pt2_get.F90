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

subroutine PT2_GET(NSIZE,LAB,VEC)

use caspt2_global, only: LUDMAT
use caspt2_module, only: CLab10, iAdr10
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSIZE
character(len=*), intent(in) :: LAB
real(kind=wp), intent(out) :: VEC(NSIZE)
integer(kind=iwp) :: I, IAD, NSZ
character(len=8) :: LAB1

I = 9-len(LAB)
if (I >= 1) then
  LAB1 = '        '
  LAB1(I:8) = LAB
else
  LAB1 = LAB(1:8)
end if

! FIND DISK ADDRESS:
do I=1,size(CLAB10)
  if (CLAB10(I) == LAB1) then
    NSZ = min(IADR10(I,2),NSIZE)
    IAD = IADR10(I,1)
    call DDAFILE(LUDMAT,2,VEC,NSZ,IAD)
#   ifdef _DEBUGPRINT_
    write(u6,*) LAB1,' SUCCESSFULLY READ FROM LUDMAT.'
    write(u6,*) '         SIZE:',NSZ,' *8 BYTES'
    write(u6,*) ' DISK ADDRESS:',IADR10(I,1)
    write(u6,'(10F12.8)') VEC(1:min(10,NSZ))
#   endif
    return
  end if
end do
write(u6,*) ' LABEL ',LAB1,' NOT FOUND ON LUDMAT.'
call ABEND()

end subroutine PT2_GET
