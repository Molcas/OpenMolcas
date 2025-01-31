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

subroutine MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDIST)
! A combination of NIGRP groups are given (IGRP)
! Find MIN and MAX for symmetry in each group
!
! Jeppe Olsen, September 1997
!              April 1998     From  MINMAX_SM_GP

use lucia_data, only: MINMAX_SM_GP

implicit none
! Input
integer NIGRP, NDIST
integer IGRP(NIGRP)
! Output
integer MNVAL(NIGRP), MXVAL(NIGRP)
! Local scratch
integer NTEST, JGRP

NTEST = 0
if (NTEST >= 100) write(6,*) ' >> Entering MINMAX_... <<'

do JGRP=1,NIGRP
  MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
  MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
end do

! Number of strings per sym and group
!do JGRP=1,NIGRP
!  call ICOPVE2(NSTSGP(1)%I,(IGRP(JGRP)-1)*NSMST+1,NSMST,LSMGP(1,JGRP))
!end do
!if (NTEST >= 1000) then
!  write(6,*) ' LSMGP'
!  call IWRTMA(LSMGP,NSMST,NIGRP,MXPOBS,NIGRP)
!end if
! Max and min sym in each group
!do JGRP=1,NIGRP
!
!  IMAX = 1
!  do ISM=1,NSMST
!    if (LSMGP(ISM,JGRP) > 0) IMAX = ISM
!  end do
!  MXVAL(JGRP) = IMAX
!
!  IMIN = NSMST
!  do ISM=NSMST,1,-1
!   if (LSMGP(ISM,JGRP) > 0) IMIN = ISM
!  end do
!  MNVAL(JGRP) = IMIN
!end do
! Total number of symmetry distributions
NDIST = 1
do JGRP=1,NIGRP
  NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
end do

if (NTEST >= 100) then
  write(6,*) ' Group combination :'
  write(6,'(5X,10I3)') (IGRP(JGRP),JGRP=1,NIGRP)
  write(6,*)
  write(6,*) ' Group Minsym Maxsym'
  write(6,*) ' ==================='
  do JGRP=1,NIGRP
    write(6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
  end do
  write(6,*)
  write(6,*) ' Total number of distributions',NDIST
end if

if (NTEST >= 1000) write(6,*) ' >> Leaving MINMAX_... <<'

end subroutine MINMAX_FOR_SYM_DIST
