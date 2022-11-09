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

subroutine WRH(LU,NSYM,NBAS,NORB,CMO,OCC,LOCC,TITLE)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LU, NSYM, NBAS(NSYM), NORB(NSYM), LOCC
real(kind=wp), intent(in) :: CMO(*), OCC(*)
character(len=*), intent(inout) :: TITLE
integer(kind=iwp) :: I, IBAS, IORB, ISYM, KCMO, KOCC, NDIV
character(len=40) :: FRMT

FRMT = '(4E20.12)'
! REWIND (LU)
KCMO = 0
NDIV = 4
if (TITLE(1:1) /= '*') TITLE = '*'//TITLE(:len(TITLE)-1)
if (locc /= 2) then
  do ISYM=1,NSYM
    do IORB=1,NORB(ISYM)
      write(LU,'(A,I5)') '* Column    ',IORB
      do IBAS=1,NBAS(ISYM),NDIV
        write(LU,FRMT) (CMO(I+KCMO),I=IBAS,min(IBAS+3,NBAS(ISYM)))
      end do
      KCMO = KCMO+NBAS(ISYM)
    end do
  end do
end if
if (LOCC == 0) return
write(LU,'(A)') Title
KOCC = 0
do ISYM=1,NSYM
  do IORB=1,NORB(ISYM),NDIV
    write(LU,FRMT) (OCC(I+KOCC),I=IORB,min(IORB+3,NORB(ISYM)))
  end do
  KOCC = KOCC+NORB(ISYM)
end do

return

end subroutine WRH
