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

subroutine TS_SYM_PNT2(IGRP,NIGRP,MAXVAL,MINVAL,ISYM,IPNT,LPNT)
! Construct pointers to start of symmetrydistributions
! for supergroup of strings with given symmetry
!
! The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
! is given as
!     1
!     +  (ISM1-MINVAL(1))
!     +  (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
!     +  (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
!     +  ...
!     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
!
! Where L is the last group of strings with nonvanishing occupation
!
! Jeppe Olsen, September 1997
!
! Version 2 : Uses IGRP and NIGRP to define supergroup

use lucia_data, only: MINMAX_SM_GP, NELFGP, NSTFSMGP
use lucia_data, only: MXPNGAS, MXPNSMST
use csm_data, only: NSMST

implicit none
integer NIGRP, ISYM, LPNT
! Specific Input
integer IGRP(NIGRP)
! Local scratch
integer ISMFGS(MXPNGAS)
!integer ITPFGS(MXPNGAS)
!-jwk-cleanup integer NELFGS(MXPNGAS)
integer NNSTSGP(MXPNSMST,MXPNGAS)
! Output
integer minval(*), maxval(*), IPNT(*)
integer NTEST, NGASL, IGAS, NBLKS, IFIRST, NSTRINT, NONEW, ISTSMM1, ISMGSN, NSTRII, IOFF, MULT, ISYMSTR

NTEST = 0
! Info on groups of strings in supergroup
NGASL = 1
do IGAS=1,NIGRP
  !ITPFGS(IGAS) = IGRP(IGAS)
  if (NELFGP(IGRP(IGAS)) > 0) NGASL = IGAS
  !. Number of strings per symmetry in each gasspace
  !call ICOPVE2(NSTSGP(1)%I,(ITPFGS(IGAS)-1)*NSMST+1,NSMST,NNSTSGP(1,IGAS))
  call ICOPVE(NSTFSMGP(1,IGRP(IGAS)),NNSTSGP(1,IGAS),NSMST)
end do

!NGASL = NIGRP
!
!do IGAS=1,NIGRP
!  do ISMST=1,NSMST
!    if (NNSTSGP(ISMST,IGAS) > 0) MAXVAL(IGAS) = ISMST
!  end do
!  do ISMST=NSMST,1,-1
!    if (NNSTSGP(ISMST,IGAS) > 0) MINVAL(IGAS) = ISMST
!  end do
!end do
do IGAS=1,NIGRP
  minval(IGAS) = MINMAX_SM_GP(1,IGRP(IGAS))
  maxval(IGAS) = MINMAX_SM_GP(2,IGRP(IGAS))
end do

if (NTEST >= 1000) then
  write(6,*) 'NIGRP:',NIGRP
  write(6,*) ' MINVAL and MAXVAL'
  call IWRTMA(MINVAL,1,NIGRP,1,NIGRP)
  call IWRTMA(MAXVAL,1,NIGRP,1,NIGRP)
end if

! Total number of symmetry blocks that will be generated
NBLKS = 1
do IGAS=1,NGASL-1
  NBLKS = NBLKS*(maxval(IGAS)-minval(IGAS)+1)
end do
if (NBLKS > LPNT) then
  write(6,*) ' Problem in TS_SYM_PNT'
  write(6,*) ' Dimension of IPNT too small'
  write(6,*) ' Actual and required length',NBLKS,LPNT
  write(6,*)
  write(6,*) ' I will Stop and wait for instructions'
  !stop ' TS_SYM_PNT too small'
  call SYSABENDMSG('lucia_util/ts_sym_pnt','Internal error','')
end if
! Loop over symmetry blocks in standard order
IFIRST = 1
NSTRINT = 0
2000 continue
if (IFIRST == 1) then
  do IGAS=1,NGASL-1
    ISMFGS(IGAS) = minval(IGAS)
  end do
else
  ! Next distribution of symmetries in NGAS -1
  call NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
  if (NONEW /= 0) goto 2001
end if
IFIRST = 0
! Symmetry of NGASL -1 spaces given, symmetry of full space
!ISTSMM1 = 1
!do IGAS=1,NGASL-1
!  call SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
!  ISTSMM1 = JSTSMM1
!end do
ISTSMM1 = ISYMSTR(ISMFGS,NGASL-1)
! sym of SPACE NGASL
call SYMCOM(2,1,ISTSMM1,ISMGSN,ISYM)
ISMFGS(NGASL) = ISMGSN
if (NTEST >= 1000) then
  write(6,*) ' next symmetry of NGASL spaces'
  call IWRTMA(ISMFGS,1,NGASL,1,NGASL)
end if
! Number of strings with this symmetry combination
NSTRII = 1
do IGAS=1,NGASL
  NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
end do
! Offset for this symmetry distribution in IOFFI
IOFF = 1
MULT = 1
do IGAS=1,NGASL-1
  IOFF = IOFF+(ISMFGS(IGAS)-minval(IGAS))*MULT
  MULT = MULT*(maxval(IGAS)-minval(IGAS)+1)
end do

IPNT(IOFF) = NSTRINT+1
NSTRINT = NSTRINT+NSTRII
if (NTEST >= 1000) write(6,*) ' IOFF, IPNT(IOFF) NSTRII ',IOFF,IPNT(IOFF),NSTRII

if (NGASL-1 > 0) goto 2000
2001 continue

if (NTEST >= 100) then
  write(6,*)
  write(6,*) ' Output from TS_SYM_PNT'
  write(6,*) ' Required total symmetry',ISYM
  write(6,*) ' Number of symmetry blocks ',NBLKS
  write(6,*)
  write(6,*) ' Offset array  for symmetry blocks'
  call IWRTMA(IPNT,1,NBLKS,1,NBLKS)
end if

end subroutine TS_SYM_PNT2
