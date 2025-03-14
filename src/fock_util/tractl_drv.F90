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
! Copyright (C) 2004, Giovanni Ghigo                                   *
!               2005, Francesco Aquilante                              *
!***********************************************************************

!***********************************************************************
! This Driver precedes the standard (TraCtl) and Cholesky (Cho_TraCtl) *
! AO/MO two-electrons transformation programs.                         *
! -------------------------------------------------------------------- *
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden                                    *
! Written:  August-September 2004                                      *
! Modified for Cholesky-MP2 May 2005                                   *
!***********************************************************************
!  TraCtl_Drv
!
!> @brief
!>   Driver that calls the Cholesky or the Conventional routine for the
!>   generation of the two-electrons integrals file in MO basis.
!> @author Giovanni Ghigo
!> @modified_by F. Aquilante
!>
!> @details
!> All programs that need the generation of the two-electrons
!> integrals file in MO basis must tell to the Cholesky routine
!> who they are through \p iType:
!>
!> - ``1``: MBPT2
!> - ``2``: CASPT2
!> - ``3``: MCLR
!> - ``4``: CC
!>
!> Programs must also tell to the Cholesky routine whether they need
!> exchange-2 integrals through logical variable \p DoExch2. Both
!> values are not used by the Conventional routine.
!>
!> @param[in] iType   Caller program
!> @param[in] DoExch2 Flag for the generation of Exch-2 integrals
!> @param[in] iPart   Partitioning of temp files
!***********************************************************************
subroutine TraCtl_Drv(iType,DoExch2,iPart)

use caspt2_global, only: CMO, LUINTM, NCMO
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iType, iPart
logical(kind=iwp), intent(in) :: DoExch2
#include "caspt2.fh"
logical(kind=iwp) :: DoCholesky
character(len=*), parameter :: SECNAM = 'TraCtl_Drv'

call DecideOnCholesky(DoCholesky)

if (DoCholesky) then

  if (iType == 1) then

    call ChoMP2_TraCtl(LUINTM,CMO,NCMO)

  ! caspt2 with cholesky does no longer use call to tractl_drv/cho_caspt2_drv
  !else if (iALGO == 0) then
  else

    call Cho_TraCtl(iType,LUINTM,CMO,NCMO,DoExch2)

  !else if (iALGO == 1) then
  !
  !  call Cho_caspt2_drv(CMO)
  !
  !else
  !
  !  call Cho_x_Quit(SecNam,' !!! Unknown algorithm !!! ',' ')

  end if

else

  call TRACTL(iPart)

end if

return

end subroutine TraCtl_Drv
