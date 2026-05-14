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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************

subroutine PadCMO(CMO,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine pads CMO's from nBas x nOrb to nBas x nBas.             *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp) :: iFrom(8), iPtr, iSym, iTo(8), ndata

!----------------------------------------------------------------------*
! Transfer orbitals.                                                   *
!----------------------------------------------------------------------*
iFrom(1) = nBas(1)*nOrb(1)
iTo(1) = nBas(1)*nOrb(1)
do iSym=1,nSym-1
  iFrom(iSym+1) = iFrom(iSym)+nBas(iSym+1)*nOrb(iSym+1)
  iTo(iSym+1) = iTo(iSym)+nBas(iSym+1)*nOrb(iSym+1)+nBas(iSym)*(nBas(iSym)-nOrb(iSym))
end do
do iSym=nSym,1,-1
  ndata = nBas(iSym)*nOrb(iSym)
  CMO(iTo(iSym)-ndata+1:iTo(iSym)) = CMO(iFrom(iSym)-ndata+1:iFrom(iSym))
  if (nBas(iSym) > nOrb(iSym)) then
    ndata = nBas(iSym)*(nBas(iSym)-nOrb(iSym))
    iPtr = iTo(iSym)+1
    CMO(iPtr:iPtr+ndata-1) = Zero
  end if
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine PadCMO
