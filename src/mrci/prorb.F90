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

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension CNO(NCMO), OCC(NBAST)
character*(LENIN8) CLEAN_BNAME
external CLEAN_BNAME

write(6,*)
call XFLUSH(6)
write(6,*) 'NATURAL ORBITALS IN AO BASIS. IN EACH SYMMETRY,'
call XFLUSH(6)
write(6,*) 'THE ORBITALS PRINTED ARE THOSE UP TO AND INCLUDING'
call XFLUSH(6)
write(6,*) 'THE LAST ORBITAL WITH OCCUPATION NUMBER LARGER'
call XFLUSH(6)
write(6,'(A,F10.7)') ' THAN THRORB = ',THRORB
call XFLUSH(6)
IEB = 0
IEM = 0
NDIV = 10
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB == 0) GO TO 100
  NPRT = 0
  do I=1,NB
    if (OCC(IEB+I) >= THRORB) NPRT = I
  end do
  if (NPRT == 0) GO TO 40
  write(6,'(/28X,''SYMMETRY LABEL'',I3)') ISYM
  call XFLUSH(6)
  do IST=1,NPRT,NDIV
    IEND = min(NPRT,IST-1+NDIV)
    write(6,'(/5X,''ORBITAL'',6X,10I8)') (I,I=IST,IEND)
    call XFLUSH(6)
    write(6,'( 5X,''OCC.NO.'',8X,10F8.5)') (OCC(IEB+I),I=IST,IEND)
    call XFLUSH(6)
    write(6,*)
    call XFLUSH(6)
    do I=1,NB
      JSMO = IEM+I+NB*(IST-1)
      JEMO = IEM+I+NB*(IEND-1)
      write(6,'(1X,I3,2X,A,10F8.4)') I,CLEAN_BNAME(NAME(IEB+I),LENIN),(CNO(J),J=JSMO,JEMO,NB)
      call XFLUSH(6)
    end do
  end do
40 continue
  IEB = IEB+NB
  IEM = IEM+NB*NB
100 continue
end do

return

end subroutine PRORB
