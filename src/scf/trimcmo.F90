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

subroutine TrimCMO(CMO1,CMO2,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine trim CMO's from nBas x nBas to nBas x nOrb.             *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CMO1(*)
real(kind=wp), intent(inout) :: CMO2(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
integer(kind=iwp) :: iFrom(8), iSym, iTo(8), ndata

!----------------------------------------------------------------------*
! Transfer orbitals.                                                   *
!----------------------------------------------------------------------*
iFrom(1) = 1
iTo(1) = 1
do iSym=1,nSym-1
  iFrom(iSym+1) = iFrom(iSym)+nBas(iSym)*nBas(iSym)
  iTo(iSym+1) = iTo(iSym)+nBas(iSym)*nOrb(iSym)
  if (iTo(iSym+1) > iFrom(iSym+1)) then
    write(u6,*) 'Error in TrimCMO'
    call Abend()
  end if
end do
do iSym=1,nSym
  ndata = nBas(iSym)*nOrb(iSym)

  ! Note that CMO1 and CMO2 might overlap. Hence, we cannot use
  ! an ordinary call to DCopy!
  ! FIXME: Overlapping non-input-only arguments are not allowed in Fortran

  if (iFrom(iSym) /= iTo(iSym)) CMO2(iTo(iSym):iTo(iSym)+nData-1) = CMO1(iFrom(iSym):iFrom(iSym)+nData-1)
end do
!----------------------------------------------------------------------*
! Finish                                                               *
!----------------------------------------------------------------------*
return

end subroutine TrimCMO
