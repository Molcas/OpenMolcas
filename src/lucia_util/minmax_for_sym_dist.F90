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
! Copyright (C) 1997,1998, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDIST)
! A combination of NIGRP groups are given (IGRP)
! Find MIN and MAX for symmetry in each group
!
! Jeppe Olsen, September 1997
!              April 1998     From  MINMAX_SM_GP

use lucia_data, only: MINMAX_SM_GP
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NIGRP, IGRP(NIGRP)
integer(kind=iwp), intent(out) :: MNVAL(NIGRP), MXVAL(NIGRP), NDIST
integer(kind=iwp) :: JGRP

#ifdef _DEBUGPRINT_
write(u6,*) ' >> Entering MINMAX_... <<'
#endif

do JGRP=1,NIGRP
  MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
  MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
end do

! Number of strings per sym and group
!do JGRP=1,NIGRP
!  LSMGP(1:NIRREP,JGRP) = NSTSGP(1)%I((IGRP(JGRP)-1)*NIRREP+1:IGRP(JGRP)*NIRREP)
!end do
!#ifdef _DEBUGPRINT_
!write(u6,*) ' LSMGP'
!call IWRTMA(LSMGP,NIRREP,NIGRP,MXPOBS,NIGRP)
!#endif
! Max and min sym in each group
!do JGRP=1,NIGRP
!
!  IMAX = 1
!  do ISM=1,NIRREP
!    if (LSMGP(ISM,JGRP) > 0) IMAX = ISM
!  end do
!  MXVAL(JGRP) = IMAX
!
!  IMIN = NIRREP
!  do ISM=NIRREP,1,-1
!    if (LSMGP(ISM,JGRP) > 0) IMIN = ISM
!  end do
!  MNVAL(JGRP) = IMIN
!end do
! Total number of symmetry distributions
NDIST = 1
do JGRP=1,NIGRP
  NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Group combination :'
write(u6,'(5X,10I3)') (IGRP(JGRP),JGRP=1,NIGRP)
write(u6,*)
write(u6,*) ' Group Minsym Maxsym'
write(u6,*) ' ==================='
do JGRP=1,NIGRP
  write(u6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
end do
write(u6,*)
write(u6,*) ' Total number of distributions',NDIST

write(u6,*) ' >> Leaving MINMAX_... <<'
#endif

end subroutine MINMAX_FOR_SYM_DIST
