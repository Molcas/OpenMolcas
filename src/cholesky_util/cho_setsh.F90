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

subroutine CHO_SETSH(IBASSH,NBASSH,NBSTSH,IBAS,NBAS,ISOSHL,NSYM,NSHELL,NBAST)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NSYM, IBAS(NSYM), NBAS(NSYM), NBAST, ISOSHL(NBAST), NSHELL
integer(kind=iwp), intent(out) :: IBASSH(NSYM,NSHELL), NBASSH(NSYM,NSHELL), NBSTSH(NSHELL)
integer(kind=iwp) :: IA, ISHL, ISYM

NBASSH(:,:) = 0
do ISYM=1,NSYM
  do IA=1,NBAS(ISYM)
    ISHL = ISOSHL(IBAS(ISYM)+IA)
    NBASSH(ISYM,ISHL) = NBASSH(ISYM,ISHL)+1
  end do
end do

do ISHL=1,NSHELL
  IBASSH(1,ISHL) = 0
  NBSTSH(ISHL) = NBASSH(1,ISHL)
  do ISYM=2,NSYM
    IBASSH(ISYM,ISHL) = NBSTSH(ISHL)
    NBSTSH(ISHL) = NBSTSH(ISHL)+NBASSH(ISYM,ISHL)
  end do
end do

end subroutine CHO_SETSH
