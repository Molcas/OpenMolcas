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

subroutine PRGASTAB(REST)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: REST(*)
integer(kind=iwp) :: IGAS, ISYM, KORB, KREST, NGAS, NSYM

write(u6,*)
write(u6,*) ' GAS restriction table printout'
write(u6,'(A,I5)') 'Table size        NSIZE=',REST(1)
write(u6,'(A,I5)') 'Table type ID     ITYPE=',REST(2)
write(u6,'(A,I5)') 'Nr of partitions  NGAS=',REST(3)
write(u6,'(A,I5)') 'Nr of symm labels NSYM =',REST(4)
write(u6,*) ' Orbital partitions:'
NGAS = REST(3)
NSYM = REST(4)
KORB = 5
write(u6,'(8X,I5,5X,8I5)') (REST(KORB+ISYM),ISYM=0,NSYM-1)
do IGAS=1,NGAS
  write(u6,'(I3,5X,I5,5X,8I5)') IGAS,(REST(KORB+ISYM+(NSYM+1)*IGAS),ISYM=0,NSYM-1)
end do
write(u6,*) ' Electron population restrictions:'
KREST = KORB+(NGAS+1)*(NSYM+1)
write(u6,'(5X,A7,5X,30I3)') 'Minimum',(REST(KREST+0+2*(IGAS-1)),IGAS=1,NGAS)
write(u6,'(5X,A7,5X,30I3)') 'Maximum',(REST(KREST+1+2*(IGAS-1)),IGAS=1,NGAS)

end subroutine PRGASTAB
