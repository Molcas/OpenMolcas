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

subroutine GenerateGeoms(Q)

use Quater_globals, only: translate, rotate, ngeoms, list, XYZ1, XYZ2
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Q(0:3)
real(kind=wp) :: Vtrans(3)
integer(kind=iwp) :: igeom, iLU
character(len=6) :: fname
integer(kind=iwp), external :: isfreeunit

do igeom=3,ngeoms+2
  list(igeom)%geo(:,:) = list(2)%geo(:,:)
end do
if (rotate) call RotateGeoms(Q)
if (translate) then
  call SetVectTrans(list(1)%nat,list(1)%geo,XYZ1,list(2)%nat,list(2)%geo,XYZ2,Vtrans)
  call TranslateGeoms(Vtrans)
end if

iLU = isfreeunit(12)
do igeom=1,ngeoms+2
  if (igeom < 100) write(fname,'(a4,i2)') 'GEOM',igeom
  if (igeom < 10) write(fname,'(a5,i1)') 'GEOM0',igeom
  call Molcas_Open(iLU,fname)
  call PrintGeom(iLU,list(igeom)%nat,list(igeom)%title,list(igeom)%geo,list(igeom)%geolbl)
  close(iLU)
end do

return

end subroutine GenerateGeoms
