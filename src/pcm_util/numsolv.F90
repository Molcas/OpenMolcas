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

function NumSolv(Solvent)

use Solvent_Data, only: Init_Solvent_Data, SolvData
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NumSolv
character(len=*), intent(in) :: Solvent
integer(kind=iwp) :: i, IdSolv
character(len=len(Solvent)) :: Solvent_

call Init_Solvent_Data()

IdSolv = 0
Solvent_ = Solvent
call Upcase(Solvent_)
do i=1,size(SolvData)
  if (Solvent_ == SolvData(i)%SName) then
    idSolv = i
    exit
  end if
end do

if (IdSolv == 0) then
  write(u6,*) ' Unrecognized solvent: ',Solvent
  write(u6,*) 'Allowed solvents are:'
  do i=1,size(SolvData)
    write(u6,*) trim(SolvData(i)%SName)
  end do
  call Abend()
end if
NumSolv = IdSolv

return

end function NumSolv
