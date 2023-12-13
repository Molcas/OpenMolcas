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

subroutine PRORB(CNO,OCC)

use mrci_global, only: BNAME, NBAS, NBAST, NCMO, NSYM, THRORB
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CNO(NCMO), OCC(NBAST)
#include "Molcas.fh"
integer(kind=iwp) :: I, IEB, IEM, IEND, IST, ISYM, J, JEMO, JSMO, NB, NDIV, NPRT
character(len=LenIn8), external :: CLEAN_BNAME

write(u6,*)
write(u6,*) 'NATURAL ORBITALS IN AO BASIS. IN EACH SYMMETRY,'
write(u6,*) 'THE ORBITALS PRINTED ARE THOSE UP TO AND INCLUDING'
write(u6,*) 'THE LAST ORBITAL WITH OCCUPATION NUMBER LARGER'
write(u6,'(A,F10.7)') ' THAN THRORB = ',THRORB
IEB = 0
IEM = 0
NDIV = 10
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB == 0) cycle
  NPRT = 0
  do I=1,NB
    if (OCC(IEB+I) >= THRORB) NPRT = I
  end do
  if (NPRT /= 0) then
    write(u6,'(/28X,''SYMMETRY LABEL'',I3)') ISYM
    do IST=1,NPRT,NDIV
      IEND = min(NPRT,IST-1+NDIV)
      write(u6,'(/5X,''ORBITAL'',6X,10I8)') (I,I=IST,IEND)
      write(u6,'( 5X,''OCC.NO.'',8X,10F8.5)') (OCC(IEB+I),I=IST,IEND)
      write(u6,*)
      do I=1,NB
        JSMO = IEM+I+NB*(IST-1)
        JEMO = IEM+I+NB*(IEND-1)
        write(u6,'(1X,I3,2X,A,10F8.4)') I,CLEAN_BNAME(BNAME(IEB+I),LenIn),(CNO(J),J=JSMO,JEMO,NB)
      end do
    end do
  end if
  IEB = IEB+NB
  IEM = IEM+NB*NB
end do

return

end subroutine PRORB
