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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

subroutine NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN,NCOMP_PER_OPEN,NCOMP)
! Number of configurations per number of open orbitals is given
! Find total number of some components, defined by number
! of components per open
!
! In practice : components are SD's, CSF's or CMB's
!
! Jeppe Olsen, Dec. 2001

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: MAXOP, NCONF_PER_OPEN(*), NCOMP_PER_OPEN(*), NCOMP
integer(kind=iwp) :: IOPEN

NCOMP = 0
do IOPEN=0,MAXOP
  NCOMP = NCOMP+NCONF_PER_OPEN(IOPEN+1)*NCOMP_PER_OPEN(IOPEN+1)
end do

end subroutine NCNF_TO_NCOMP
