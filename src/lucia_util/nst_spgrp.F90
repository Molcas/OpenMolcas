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
! Copyright (C) 1997, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NST_SPGRP(NGRP,IGRP,ISM_TOT,NSTSGP,NSMST,NSTRIN,NDIST)
! Number of strings for given combination of groups and symmetry.
!
! Input
!
!   NGRP : Number of active groups
!   IGRP : The active groups
!   ISM_TOT : Total symmetry of supergroup
!   NSTSGP  : Number of strings per symmetry and supergroup
!   NSMST   : Number of string symmetries
!
! Output
!
!   NSTRIN : Number of strings with symmetry ISM_TOT
!   NDIST  : Number of symmetry distributions
!
! Jeppe Olsen, September 1997

use lucia_data, only: INGRP_VAL, ISMDFGP, ISMSCR, MXPNGAS, NACTSYM
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NGRP, IGRP(NGRP), ISM_TOT, NSMST, NSTSGP(NSMST,*)
integer(kind=iwp), intent(out) :: NSTRIN, NDIST
integer(kind=iwp) :: IFIRST, ISM(MXPNGAS), JGRP, LDIST, LENGTH, MNSM(MXPNGAS), MXSM(MXPNGAS), NDISTX, NONEW
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I

write(u6,*) ' ===================='
write(u6,*) ' NST_SPGP is speaking'
write(u6,*) ' ===================='

write(u6,*) ' Supergroup in action :'
write(u6,'(A,I3  )') ' Number of active spaces ',NGRP
write(u6,'(A,20I3)') ' The active groups       ',(IGRP(I),I=1,NGRP)
#endif
! Set up min and max values for symmetries
call MINMAX_FOR_SYM_DIST(NGRP,IGRP,MNSM,MXSM,NDISTX)
! Loop over symmetry distributions
IFIRST = 1
LENGTH = 0
NDIST = 0
do
  ! Next symmetry distribution
  !GLM call NEXT_SYM_DISTR(NGRP,MNSM,MXSM,ISM,ISM_TOT,IFIRST,NONEW)
  !call NEXT_SYM_DISTR(NGASL,MNVLK,MXVLK,ISMFGS,KSM,KFIRST,NONEW)
  ! GLMJ Giovanni Li Manni modification  Feb/March 2012
  call NEXT_SYM_DISTR_NEW(NSMST,INGRP_VAL,IGRP,NGRP,ISM,ISM_TOT,IFIRST,NONEW,ISMDFGP,NACTSYM,ISMSCR)
  if (NONEW /= 0) exit
  LDIST = 1
  do JGRP=1,NGRP
    LDIST = LDIST*NSTSGP(ISM(JGRP),IGRP(JGRP))
  end do
  LENGTH = LENGTH+LDIST
  NDIST = NDIST+1
end do

NSTRIN = LENGTH

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of strings obtained ',LENGTH
write(u6,*) ' Number of symmetry-distributions',NDIST
#endif

end subroutine NST_SPGRP
