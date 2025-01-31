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

subroutine SMDFGP_GEN(NGRP,NSMST,MXPNS,NSTFSMGP,NACTSYM,ISMDFGP)

implicit real*8(A-H,O-Z)
integer ISMDFGP(NSMST,NGRP)
integer NSTFSMGP(MXPNS,NGRP)
integer NACTSYM(NGRP)
integer IOFF

NTEST = 0

if (NTEST >= 1000) then
  write(6,*) 'Entering  SMDFGP_GEN'
  write(6,*) 'NGRP : ',NGRP
  write(6,*) 'NSMST: ',NSMST
  write(6,*) 'MXPNS',MXPNS
  write(6,*) 'INPUT: NSTFSMGP'
  do ISYM=1,NSMST
    write(6,'(40I5)') (NSTFSMGP(ISYM,IGRP),IGRP=1,NGRP)
  end do
end if

do IGRP=1,NGRP
  IOFF = 0
  NACTSYM(IGRP) = IOFF
  do ISYM=1,NSMST
    ISMDFGP(ISYM,IGRP) = 0
    if (NSTFSMGP(ISYM,IGRP) /= 0) then
      IOFF = IOFF+1
      ISMDFGP(IOFF,IGRP) = ISYM
    end if
  end do
  NACTSYM(IGRP) = IOFF
end do

if (NTEST >= 100) then
  write(6,*) 'Number of Active Symm per GRP:'
  write(6,*) (NACTSYM(IGRP),IGRP=1,NGRP)
  write(6,*) 'Symmetries allowed by each group:'
  do ISYM=1,NSMST
    write(6,'(40I2)') (ISMDFGP(ISYM,IGRP),IGRP=1,NGRP)
  end do
end if

end subroutine SMDFGP_GEN
