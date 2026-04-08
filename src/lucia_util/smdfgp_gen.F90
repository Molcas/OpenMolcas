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

!#define _DEBUGPRINT_
subroutine SMDFGP_GEN(NGRP,NSMST,MXPNS,NSTFSMGP,NACTSYM,ISMDFGP)

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NGRP, NSMST, MXPNS, NSTFSMGP(MXPNS,NGRP)
integer(kind=iwp), intent(out) :: NACTSYM(NGRP), ISMDFGP(NSMST,NGRP)
integer(kind=iwp) :: IGRP, IOFF, ISYM

#ifdef _DEBUGPRINT_
write(u6,*) 'Entering  SMDFGP_GEN'
write(u6,*) 'NGRP : ',NGRP
write(u6,*) 'NSMST: ',NSMST
write(u6,*) 'MXPNS',MXPNS
write(u6,*) 'INPUT: NSTFSMGP'
do ISYM=1,NSMST
  write(u6,'(40I5)') (NSTFSMGP(ISYM,IGRP),IGRP=1,NGRP)
end do
#endif

do IGRP=1,NGRP
  IOFF = 0
  do ISYM=1,NSMST
    ISMDFGP(ISYM,IGRP) = 0
    if (NSTFSMGP(ISYM,IGRP) /= 0) then
      IOFF = IOFF+1
      ISMDFGP(IOFF,IGRP) = ISYM
    end if
  end do
  NACTSYM(IGRP) = IOFF
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'Number of Active Symm per GRP:'
write(u6,*) (NACTSYM(IGRP),IGRP=1,NGRP)
write(u6,*) 'Symmetries allowed by each group:'
do ISYM=1,NSMST
  write(u6,'(40I2)') (ISMDFGP(ISYM,IGRP),IGRP=1,NGRP)
end do
#endif

end subroutine SMDFGP_GEN
