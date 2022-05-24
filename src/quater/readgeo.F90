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
!  readgeo
!
!> @brief
!>   Reads a geometry and stores it into memory
!> @author Y. Carissan
!>
!> @details
!> Reads a geometry in the XYZ format and stores it
!> in memory at \c list(ig)%geo.
!>
!> @param[in] iLU logic unit number
!> @param[in] ig  geometry index
!***********************************************************************

subroutine readgeo(iLU,ig)

use Quater_globals, only: debug, list
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iLU, ig
integer(kind=iwp) :: natoms, iat
character(len=180) :: Line
character(len=20) :: lbl
character(len=6) :: cName
character(len=180), external :: Get_ln

if ((ig < 1).or.(ig > 2)) call SysAbendMsg('ReadGeo','Wrong ig ','Shoot the programmer')

! read natoms
Line = Get_ln(iLU)
call Put_ln(Line)
call Get_I1(1,natoms)

if (natoms > 500) call SysAbendMsg('ReadGeo','Too many atoms in geom','')

if (debug) write(u6,*) 'In READGEO : Nat=',natoms

list(ig)%nat = natoms

! allocate memory for the geometry

if (ig < 10) then
  write(cName,'(a5,i1)') 'GEOM0',ig
  list(ig)%title = cName
else if (ig < 100) then
  write(cName,'(a4,i2)') 'GEOM',ig
  list(ig)%title = cName
end if
call mma_allocate(list(ig)%geo,3,list(ig)%nat,label=cName)
call mma_allocate(list(ig)%geolbl,list(ig)%nat,label=cName)

! read title

Line = Get_ln(iLU)
list(ig)%title = trim(line)

! read label and coords

iat = 0
do
  Line = Get_ln(iLU)
  call Put_ln(Line)
  call Get_S(1,lbl,1)
  if (lbl == 'END ') then
    exit
  else
    iat = iat+1
    if (iat > list(ig)%nat) call SysAbendMsg('ReadGeo','More atoms read than declared','')
    list(ig)%geolbl(iat) = lbl
    call Get_F(2,list(ig)%geo(:,iat),3)
  end if
end do

call PrintGeom(-1_iwp,list(ig)%nat,list(ig)%title,list(ig)%geo,list(ig)%geolbl)

return

end subroutine readgeo
