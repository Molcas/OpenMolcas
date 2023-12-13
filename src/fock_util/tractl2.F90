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

subroutine TraCtl2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)
!***********************************************************************
!                                                                      *
!     main control section for Cholesky-based vs conventional          *
!     - transformation of ERIs from AO to MO basis                     *
!     - Fock matrix generation                                         *
!                                                                      *
!***********************************************************************

use Fock_util_global, only: ALGO, DoCholesky
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*), D1I(*), D1A(*), ExFac
real(kind=wp), intent(inout) :: PUVX(*), FI(*), FA(*)
real(kind=wp), intent(_OUT_) :: TUVX(*)
integer(kind=iwp), intent(in) :: IPR
logical(kind=iwp), intent(in) :: lSquare
#include "rasdim.fh"
#include "general.fh"
#include "wadr.fh"
integer(kind=iwp) :: iDisk, irc
logical(kind=iwp) :: TraOnly

!call DecideOnCholesky(DoCholesky)

!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

if (.not. DoCholesky) then

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  call TRA_CTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

else if (ALGO == 1) then

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  TraOnly = .false.
  call CHO_CAS_DRV(irc,CMO,D1I,FI,D1A,FA,PUVX,TraOnly)

# ifdef _MOLCAS_MPP_
  ! --------------------------------------------------
  ! Synchronize Fock matrices if running parallel:
  if ((nProcs > 1) .and. Is_Real_Par()) then
    call GADsum(FI,nTot1)
    call GADsum(FA,nTot1)
    ! Synchronize PUVX if running parallel:
    call GADsum(PUVX,nPWXY)
  end if
# endif
  ! --------------------------------------------------
  ! select integrals TUVX
  call Get_TUVX(PUVX,TUVX)
  ! save integrals on disk
  ! nPWXY is computed in cho_eval_waxy and stored in wadr.fh
  iDisk = 0
  call DDaFile(LUINTM,1,PUVX,nPWXY,iDisk)

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

else if (ALGO == 2) then

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  TraOnly = .false.
  call CHO_CAS_DRV(irc,CMO,D1I,FI,D1A,FA,PUVX,TraOnly)

  if (irc /= 0) then
    write(u6,*) 'TRACTL2: Cho_cas_drv non-Zero return code. rc= ',irc
    call Abend()
  end if

  ! Synchronization for parallel runs is done in cho_cas_drv

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

end if

!)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

return

end subroutine TraCtl2
