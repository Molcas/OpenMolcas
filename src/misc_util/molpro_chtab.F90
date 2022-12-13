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
! Copyright (C) 2014, Naoki Nakatani                                   *
!***********************************************************************

subroutine Molpro_ChTab(nIrrep,Label,iChMolpro)
!***********************************************************************
!                                                                      *
! Object: To convert MOLCAS character table to MOLPRO format           *
!         This is a part of integral_util/chtab.f                      *
!                                                                      *
!         Written by N.Nakatani Nov. 2014                              *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nIrrep
character(len=3), intent(out) :: Label
integer(kind=iwp), intent(out) :: iChMolpro(8)
integer(kind=iwp) :: i, iOper(8)
logical(kind=iwp) :: Rot

!***********************************************************************
call Get_iArray('Symmetry operations',iOper,nIrrep)

iChMolpro(:) = 0

if (nIrrep == 1) then
  !***** C1  symmetry ******
  Label = 'c1 '
  iChMolpro(1) = 1
else if (nIrrep == 2) then
  if (iOper(2) == 7) then
    !***** Ci  symmetry ******
    Label = 'ci '
  else if ((iOper(2) == 1) .or. (iOper(2) == 2) .or. (iOper(2) == 4)) then
    !***** Cs  symmetry ******
    Label = 'cs '
  else
    !***** C2  symmetry ******
    Label = 'c2 '
  end if
  iChMolpro(1) = 1
  iChMolpro(2) = 2
else if (nIrrep == 4) then
  if ((iOper(2) == 7) .or. (iOper(3) == 7) .or. (iOper(4) == 7)) then
    !***** C2h symmetry ******
    Label = 'c2h'
    ! MOLCAS::[ ag, bg, au, bu ] => MOLPRO::[ ag, au, bu, bg ]
    iChMolpro(1) = 1
    iChMolpro(2) = 4
    iChMolpro(3) = 2
    iChMolpro(4) = 3
  else
    Rot = .true.
    do i=1,nIrrep
      if ((iOper(i) == 1) .or. (iOper(i) == 2) .or. (iOper(i) == 4)) Rot = .false.
    end do
    if (Rot) then
      !***** D2  symmetry ******
      Label = 'd2 '
      ! MOLCAS::[ a, b2, b1, b3 ] => MOLPRO::[ a, b3, b2, b1 ]
      iChMolpro(1) = 1
      iChMolpro(2) = 3
      iChMolpro(3) = 4
      iChMolpro(4) = 2
    else
      !***** C2v symmetry ******
      Label = 'c2v'
      ! MOLCAS::[ a1, b1, a2, b2 ] => MOLPRO::[ a1, b1, b2, a2 ]
      iChMolpro(1) = 1
      iChMolpro(2) = 2
      iChMolpro(3) = 4
      iChMolpro(4) = 3
    end if
  end if
else if (nIrrep == 8) then
  !***** D2h symmetry ******
  Label = 'd2h'
  ! MOLCAS::[ ag, b2g, b1g, b3g, au, b2u, b1u, b3u ] => MOLPRO::[ ag, b3u, b2u, b1g, b1u, b2g, b3g, au ]
  iChMolpro(1) = 1
  iChMolpro(2) = 6
  iChMolpro(3) = 4
  iChMolpro(4) = 7
  iChMolpro(5) = 8
  iChMolpro(6) = 3
  iChMolpro(7) = 5
  iChMolpro(8) = 2
else
  call WarningMessage(2,'MOLPRO_ChTab: Illegal value of nIrrep')
  write(u6,*) 'nIrrep=',nIrrep
  call Abend()
end if

return

end subroutine Molpro_ChTab
