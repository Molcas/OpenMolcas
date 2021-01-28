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
!  RdInput_Quater
!
!> @brief
!>   Reads the input of the quater program
!> @author Y. Carissan
!>
!> @param[out] U1 vector needed for the rotation
!> @param[out] U2 vector needed for the rotation
!> @param[out] V1 vector needed for the rotation
!> @param[out] V2 vector needed for the rotation
!>
!> @details
!> Reads the input of the quater program.
!***********************************************************************

subroutine RdInput_Quater(U1,U2,V1,V2)

use Quater_globals, only: debug, rotate, translate, ngeoms, list, XYZ1, XYZ2
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: U1(3), U2(3), V1(3), V2(3)
integer(kind=iwp) :: iLU, ig
character(len=6) :: cName
character(len=180) :: Line
character(len=20) :: key
logical(kind=iwp) :: AxisSet, NewAxisSet, Geo1Set, Geo2Set, XYZ1Set, XYZ2Set
character(len=180), external :: Get_ln

AxisSet = .false.
NewAxisSet = .false.
Geo1Set = .false.
Geo2Set = .false.
XYZ1Set = .false.
XYZ2Set = .false.
rotate = .true.
translate = .true.
ngeoms = 1
call SpoolInp(iLU)

! skip &quater &end

Line = Get_ln(iLU)
do
  Line = Get_ln(iLU)
  call Put_ln(Line)
  call Get_S(1,key,1)
  if (debug) write(u6,*) 'KEY:',key
  if (key(1:4) == 'AXIS') then
    Line = Get_ln(iLU)
    call Put_ln(Line)
    call Get_F(1,U1,3)
    Line = Get_ln(iLU)
    call Put_ln(Line)
    call Get_F(1,U2,3)
    AxisSet = .true.
  else if (key(1:4) == 'NEWA') then
    Line = Get_ln(iLU)
    call Put_ln(Line)
    call Get_F(1,V1,3)
    Line = Get_ln(iLU)
    call Put_ln(Line)
    call Get_F(1,V2,3)
    NewAxisSet = .true.
  else if (key(1:4) == 'DEBU') then
    debug = .true.
  else if (key(1:4) == 'GEO1') then
    call readgeo(iLU,1)
    Geo1Set = .true.
  else if (key(1:4) == 'GEO2') then
    call readgeo(iLU,2)
    Geo2Set = .true.
  else if (key(1:4) == 'XYZ1') then
    Line = Get_ln(iLU)
    call Put_ln(Line)
    call Get_I(1,XYZ1,3)
    XYZ1Set = .true.
  else if (key(1:4) == 'XYZ2') then
    Line = Get_ln(iLU)
    call Put_ln(Line)
    call Get_I(1,XYZ2,3)
    XYZ2Set = .true.
  !else if (key(1:4) == 'NGEO') then
  !  Line = Get_ln(iLU)
  !  call Put_ln(Line)
  !  call Get_I(1,ngeoms,1)
  else if (key(1:4) == 'NOTR') then
    translate = .false.
  else if (key(1:4) == 'NORO') then
    rotate = .false.
  else if (key(1:4) == 'END ') then
    exit
  else
    call SysAbendMsg('RdInput_Quater','Keyword not relevant : ',key)
  endif
end do

if (.not.AxisSet) then
  if (.not.XYZ1Set) then
    call SysAbendMsg('RdInput_Quater','Reference Axis not set','AXIS or XYZ1 Keyword mandatory')
  end if
else
  if (XYZ1Set) then
    call SysAbendMsg('RdInput_Quater','Reference Axis not set properly','AXIS and XYZ1 Keywords are exclusive')
  end if
end if

if (.not.NewAxisSet) then
  if (.not.XYZ2Set) then
    call SysAbendMsg('RdInput_Quater','New Axis not set :','NEWAXIS or XYZ2 Keyword mandatory')
  end if
else
  if (XYZ2Set) then
    call SysAbendMsg('RdInput_Quater','New Axis not set properly','NEWAXIS and XYZ2 Keywords are exclusive')
  end if
end if

if (XYZ1Set.and..not.GEO1Set) then
  call SysAbendMsg('RdInput_Quater','XYZ1 keyword requires GEO1 definition','')
end if

if (XYZ2Set.and..not.GEO2Set) then
  call SysAbendMsg('RdInput_Quater','XYZ2 keyword requires GEO2 definition','')
end if

if (translate.and..not.(XYZ1Set.and.XYZ2Set)) then
  call SysAbendMsg('RdInput_Quater','Translation cannot be done if both','XYZ1 and XYZ2 are not set')
end if

if (XYZ1Set) call SetVect(list(1)%nat,list(1)%geo,XYZ1,V1,V2)
if (XYZ2Set) call SetVect(list(2)%nat,list(2)%geo,XYZ2,U1,U2)

do ig=3,ngeoms+2
  if (ig < 10) then
    write(cName,'(a5,i1)') 'GEOM0',ig
    list(ig)%title = cName
  else if (ig < 100) then
    write(cName,'(a4,i2)') 'GEOM',ig
    list(ig)%title = cName
  end if
  list(ig)%nat = list(2)%nat
  call mma_allocate(list(ig)%geo,3,list(ig)%nat,label=cName)
  call mma_allocate(list(ig)%geolbl,list(ig)%nat,label=cName)
  list(ig)%geolbl(:) = list(2)%geolbl(:)
end do

return

end subroutine RdInput_Quater
