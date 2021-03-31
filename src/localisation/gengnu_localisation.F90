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

subroutine GenGnu_Localisation(FilNam,Den,Coord,n)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
character(len=12), intent(in) :: FilNam
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: Den(n,n), Coord(3,n)
integer(kind=iwp) :: i, j, Lu
real(kind=wp) :: Dist, D, Fac
character(len=16) :: FName
character(len=4) :: Postfix = '.dat'
integer(kind=iwp), external :: isFreeUnit

! Open data file.
! ---------------

Lu = isFreeUnit(11)
write(FName,'(A12,A4)') FilNam,Postfix
if (FName(1:1) == ' ') FName(1:1) = 'G'
do i=2,12
  if (FName(i:i) == ' ') FName(i:i) = '_'
end do
call Molcas_Open(Lu,FName)
rewind(Lu)

! Write data file.
! ----------------

Fac = One/log(10.0_wp)
do j=1,n
  do i=j,n
    Dist = sqrt((Coord(1,j)-Coord(1,i))**2+(Coord(2,j)-Coord(2,i))**2+(Coord(3,j)-Coord(3,i))**2)
    if (abs(Den(i,j)) < 1.0e-16_wp) then
      D = -40.0_wp
    else
      D = Fac*log(abs(Den(i,j)))
    end if
    write(Lu,'(1X,1P,D20.10,1X,D20.10)') Dist,D
  end do
end do

! Close data file.
! ----------------

close(Lu,status='Keep')

end subroutine GenGnu_Localisation
