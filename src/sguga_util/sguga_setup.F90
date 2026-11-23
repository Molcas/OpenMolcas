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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2026, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
! Initialization of a SGUGA state (SG_Init, SG_Mk_SDREO).
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_Setup

use SGUGA_Core
use SGUGA_Tables
use SGUGA_Dets

#include "intent.fh"

implicit none

contains

subroutine SG_Init(iState,nSym,nActEl,iSpin, &
                   nRas,nRasEl,nRsPrt, &
                   xLevel,xL2Act,xNLEV,xNSM,xnRasMx)

  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xnLev, xNSM(MxLev), xnRasMx(nRsPrt)

  type(TRStruct) :: TRS

  call SG_Init_Light(iState,nSym,nActEl,iSpin, &
                     nRas,nRasEl,nRsPrt, &
                     xLevel,xL2Act,xnLev,xNSM,xnRasMx=xnRasMx)

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

  call MKMAW(SGS(iState))

! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.

  call MKSEG(SGS(iState),CIS(iState),EXS(iState))

  !Create the transition infrastructure.
  call MkTrans(SGS(iState),CIS(iState),TRS)

  ! Count coupling coefficients in compressed blocks indexed by excitation/operator type, symmetry and midvertex.
  call MkNRCOUP(SGS(iState),CIS(iState),EXS(iState),TRS)

  ! Explicitly generates the compressed coupling tuples '(left walk, right walk, value index)' and compacts repeated
  ! numerical values into 'VTab'.
  call MKCOUP(SGS(iState),CIS(iState),EXS(iState),TRS)

  call Trans_Free(TRS)

  ! Generate flat tables for sg_Epq_Psi_flat

  call MK_FLAT_Epq(SGS(iState))

  ! Generate the Configuration-to-Determinant Coefficients
  ! Generate Occupation Configurations a la Lucia.
  ! Temporarily just enabled for CAS
  ! (CDCInfo is needed by SG_ReOrd. The determinant reordering SDREO is
  ! only needed by the SD/CSF codec and is built on first use by
  ! SG_Mk_SDREO: building it here made e.g. a CAS(12,12) RASSCF about 75%
  ! slower.)

  if (SGS(iState)%nRsPrt == 1) then
    call Mk_CDC(SGS(iState))
  end if

end subroutine SG_Init

!> Build the determinant reordering (SDREO) for state iState, if not yet
!> done. It is needed only by the SD/CSF codec (CASVB with Faroald, and the
!> MOLCAS_CI_VERIFY=SIGMA checks). Only available for CAS, where SG_Init has
!> already built the CDCInfo tables it relies on.
subroutine SG_Mk_SDREO(iState)
  integer(kind=iwp), intent(in) :: iState
  if (allocated(SGS(iState)%SDREO)) return
  if (SGS(iState)%nRsPrt /= 1) then
    write(u6,*) 'SG_Mk_SDREO: the SD/CSF codec is only available for CASSCF'
    call Abend()
  end if
  if (.not. allocated(SGS(iState)%CDCInfo)) call Mk_CDC(SGS(iState))
  call Mk_SDREO_SGUGA(SGS(iState))
end subroutine SG_Mk_SDREO

end module SGUGA_Setup
