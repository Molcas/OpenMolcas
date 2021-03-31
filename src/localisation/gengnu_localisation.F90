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

implicit none
character*12 FilNam
integer n
real*8 Den(n,n), Coord(3,n)

integer isFreeUnit
external isFreeUnit

character*4 Postfix
character*16 FName
parameter(Postfix='.dat')

integer i, j, Lu
real*8 Dist, D, Fac

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

Fac = 1.0d0/log(1.0d1)
do j=1,n
  do i=j,n
    Dist = sqrt((Coord(1,j)-Coord(1,i))**2+(Coord(2,j)-Coord(2,i))**2+(Coord(3,j)-Coord(3,i))**2)
    if (abs(Den(i,j)) < 1.0d-16) then
      D = -4.0d1
    else
      D = Fac*log(abs(Den(i,j)))
    end if
    write(Lu,'(1X,1P,D20.10,1X,D20.10)') Dist,D
  end do
end do

! Close data file.
! ----------------

close(Lu,Status='Keep')

end subroutine GenGnu_Localisation
