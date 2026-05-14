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
subroutine TS_SYM_PNT2(IGRP,NIGRP,MXVAL,MNVAL,ISYM,IPNT,LPNT)
! Construct pointers to start of symmetry distributions
! for supergroup of strings with given symmetry
!
! The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
! is given as
!     1
!     +  (ISM1-MNVAL(1))
!     +  (ISM2-MNVAL(2))*(MXVAL(1)-MNVAL(1)+1)
!     +  (ISM3-MNVAL(3))*(MXVAL(1)-MNVAL(1)+1)*(MXVAL(2)-MNVAL(2)+1)
!     +  ...
!     +  (ISM L-1-MNVAL(L-1))*Prod(i=1,L-2)(MXVAL(i)-MNVAL(i)+1)
!
! Where L is the last group of strings with nonvanishing occupation
!
! Jeppe Olsen, September 1997
!
! Version 2 : Uses IGRP and NIGRP to define supergroup

use Symmetry_Info, only: Mul
use lucia_data, only: MINMAX_SM_GP, MXPNGAS, MXPNSMST, NELFGP, NIRREP, NSTFSMGP
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NIGRP, IGRP(NIGRP), ISYM, LPNT
integer(kind=iwp), intent(out) :: MNVAL(NIGRP), MXVAL(NIGRP), IPNT(LPNT)
integer(kind=iwp) :: IFIRST, IGAS, IOFF, ISMFGS(MXPNGAS), ISTSMM1, ISYMSTR, MULT, NBLKS, NGASL, NNSTSGP(MXPNSMST,MXPNGAS), NONEW, &
                     NSTRII, NSTRINT

! Info on groups of strings in supergroup
NGASL = 1
do IGAS=1,NIGRP
  !ITPFGS(IGAS) = IGRP(IGAS)
  if (NELFGP(IGRP(IGAS)) > 0) NGASL = IGAS
  ! Number of strings per symmetry in each gasspace
  !NNSTSGP(1:NIRREP,IGAS) = NSTSGP(1)%I((ITPFGS(IGAS)-1)*NIRREP+1:ITPFGS(IGAS)*NIRREP)
  NNSTSGP(1:NIRREP,IGAS) = NSTFSMGP(1:NIRREP,IGRP(IGAS))
end do

!NGASL = NIGRP

!do IGAS=1,NIGRP
!  do ISMST=1,NIRREP
!    if (NNSTSGP(ISMST,IGAS) > 0) MXVAL(IGAS) = ISMST
!  end do
!  do ISMST=NIRREP,1,-1
!    if (NNSTSGP(ISMST,IGAS) > 0) MNVAL(IGAS) = ISMST
!  end do
!end do
do IGAS=1,NIGRP
  MNVAL(IGAS) = MINMAX_SM_GP(1,IGRP(IGAS))
  MXVAL(IGAS) = MINMAX_SM_GP(2,IGRP(IGAS))
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'NIGRP:',NIGRP
write(u6,*) ' MNVAL and MXVAL'
call IWRTMA(MNVAL,1,NIGRP,1,NIGRP)
call IWRTMA(MXVAL,1,NIGRP,1,NIGRP)
#endif

! Total number of symmetry blocks that will be generated
NBLKS = product(MXVAL(1:NGASL-1)-MNVAL(1:NGASL-1)+1)
if (NBLKS > LPNT) then
  write(u6,*) ' Problem in TS_SYM_PNT'
  write(u6,*) ' Dimension of IPNT too small'
  write(u6,*) ' Actual and required length',NBLKS,LPNT
  write(u6,*)
  write(u6,*) ' I will Stop and wait for instructions'
  !stop ' TS_SYM_PNT too small'
  call SYSABENDMSG('lucia_util/ts_sym_pnt','Internal error','')
end if
! Loop over symmetry blocks in standard order
IFIRST = 1
NSTRINT = 0
do
  if (IFIRST == 1) then
    ISMFGS(1:NGASL-1) = MNVAL(1:NGASL-1)
  else
    ! Next distribution of symmetries in NGAS -1
    call NXTNUM3(ISMFGS,NGASL-1,MNVAL,MXVAL,NONEW)
    if (NONEW /= 0) exit
  end if
  IFIRST = 0
  ! Symmetry of NGASL -1 spaces given, symmetry of full space
  !ISTSMM1 = 1
  !do IGAS=1,NGASL-1
  !  ISTSMM1 = Mul(ISTSMM1,ISMFGS(IGAS))
  !end do
  ISTSMM1 = ISYMSTR(ISMFGS,NGASL-1)
  ! sym of SPACE NGASL
  ISMFGS(NGASL) = Mul(ISTSMM1,ISYM)
# ifdef _DEBUGPRINT_
  write(u6,*) ' next symmetry of NGASL spaces'
  call IWRTMA(ISMFGS,1,NGASL,1,NGASL)
# endif
  ! Number of strings with this symmetry combination
  NSTRII = 1
  do IGAS=1,NGASL
    NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
  end do
  ! Offset for this symmetry distribution in IOFFI
  IOFF = 1
  MULT = 1
  do IGAS=1,NGASL-1
    IOFF = IOFF+(ISMFGS(IGAS)-MNVAL(IGAS))*MULT
    MULT = MULT*(MXVAL(IGAS)-MNVAL(IGAS)+1)
  end do

  IPNT(IOFF) = NSTRINT+1
  NSTRINT = NSTRINT+NSTRII
# ifdef _DEBUGPRINT_
  write(u6,*) ' IOFF, IPNT(IOFF) NSTRII ',IOFF,IPNT(IOFF),NSTRII
# endif

  if (NGASL <= 1) exit
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' Output from TS_SYM_PNT'
write(u6,*) ' Required total symmetry',ISYM
write(u6,*) ' Number of symmetry blocks ',NBLKS
write(u6,*)
write(u6,*) ' Offset array  for symmetry blocks'
call IWRTMA(IPNT,1,NBLKS,1,NBLKS)
#endif

end subroutine TS_SYM_PNT2
