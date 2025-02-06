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

subroutine TODSCN(VEC,NREC,LREC,LBLK,LU)
! Write VEC as multiple record file accordin to NREC and LREC

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: VEC(*)
integer(kind=iwp) :: NREC, LREC(NREC), LBLK, LU
integer(kind=iwp) :: IOFF, IREC

IOFF = 1
do IREC=1,NREC
  !write(u6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
  !write(u6,*) ' Input record'
  !call WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
  if (LREC(IREC) >= 0) then
    call ITODS(LREC(IREC),1,LBLK,LU)
    call TODSC(VEC(IOFF),LREC(IREC),LBLK,LU)
    IOFF = IOFF+LREC(IREC)
  else
    call ITODS(-LREC(IREC),1,LBLK,LU)
    call ZERORC(LU,0)
  end if
end do

end subroutine TODSCN
