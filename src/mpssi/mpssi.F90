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
! Copyright (C) 2018, Stefan Knecht                                    *
!***********************************************************************

subroutine mpssi(iReturn)

#ifdef _DMRG_
use qcmaquis_interface_cfg, only: doDMRG
#endif
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
! ----------------------------------------------------------------------

!> set QCMaquis as default (actually the only possible) DMRG driver in RASSI
#ifdef _DMRG_
doDMRG = .true.
#endif

!> print info about MPS-SI
write(u6,'(/a/,a//,a/,a/,a/,a/,a/)')                                  &
      '   ---------------------------------------------------------', &
      '   Matrix-Product-State State-Interaction (MPS-SI) activated', &
      '   Please cite for the MPS-SI framework:',                     &
      '   S. Knecht, S. Keller, J. Autschbach, M. Reiher,',           &
      '   J. Chem. Theory Comput., 12, 5881-5894 (2016).',            &
      '   ---------------------------------------------------------'

!> call state-interaction workhorse aka RASSI
iReturn = 0
call rassi(iReturn)

#ifdef _DMRG_
!> reset in case we call RASSI afterwards requesting a CI driver
doDMRG = .false.
#endif

end subroutine mpssi
