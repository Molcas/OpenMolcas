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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************

subroutine ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,printMOext)
! Author: Y. Carissan

use Localisation_globals, only: Debug
use Constants, only: Zero, One, Two, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: Functional
logical(kind=iwp), intent(in) :: printMOext
integer(kind=iwp) :: iAt, iMO_s
real(kind=wp) :: d_s
character(len=30) :: MOtype
real(kind=wp), parameter :: thr = 0.1_wp

if (Debug) then
  write(u6,*) 'In ComputeFunc:'
  write(u6,*) '-----------------'
end if

MOtype = 'delocalized'

Functional = Zero
do iMO_s=1,nOrb2Loc
  ! calculating the Pipek measure of localization, that tells over how many atoms, MO s extends,
  d_s = Zero
  do iAt=1,nAtoms
    d_s = d_s+PA(iMO_s,iMO_s,iAt)**2
  end do
  Functional = Functional+d_s

  if (d_s > 1.0e-12_wp) then
    d_s = One/d_s
  else
    if (printMOext) write(u6,*) 'WARNING: d_s^{-1} is zero for MO ',iMO_s
    d_s = huge(d_s)
  end if

  if (d_s < One+thr) then
    MOtype = 'core orbital or lone pair'
  else if (abs(Two-d_s) < thr) then
    MOtype = 'bond orbital'
  else if (d_s > OneHalf+thr) then
    MOtype = 'delocalized'
  end if

  if (printMOext) write(u6,'(A,I4,A,F8.3,1X,A,A)') ' MO',iMO_s,')   ',d_s,' atoms  -> ',MOtype
end do

end subroutine ComputeFunc
