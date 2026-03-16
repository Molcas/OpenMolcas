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

subroutine PRPRTTAB(IPRTTAB)

use Definitions, only: u6

implicit none
integer IPRTTAB(*)
integer IPART, NPART
integer NSYM, ISYM

write(u6,*)
write(u6,*) ' Partition table printout'
write(u6,'(A,I5)') 'Table size        NSIZE=',IPRTTAB(1)
write(u6,'(A,I5)') 'Table type ID     ITYPE=',IPRTTAB(2)
write(u6,'(A,I5)') 'Nr of partitions  NPART=',IPRTTAB(3)
write(u6,'(A,I5)') 'Nr of symm labels NSYM =',IPRTTAB(4)
NPART = IPRTTAB(3)
NSYM = IPRTTAB(4)
write(u6,'(8X,I5,5X,8I5)') (IPRTTAB(u6+ISYM),ISYM=0,NSYM)
do IPART=1,NPART
  write(u6,'(I3,5X,I5,5X,8I5)') IPART,(IPRTTAB(5+ISYM+(NSYM+1)*IPART),ISYM=0,NSYM)
end do

end subroutine PRPRTTAB
