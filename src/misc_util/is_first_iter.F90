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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

logical function Is_First_Iter()

implicit none
character(Len=80) :: EnvStr
integer :: Iter, Iter_S
logical :: Found
integer, dimension(7) :: Tmp

! If this is the first iteration in a "saddle" branch
call qpg_iScalar('Saddle Iter',Found)
if (Found) then
  call Get_iScalar('Saddle Iter',Iter_S)
  if (Iter_S == 0) then
    Is_First_Iter = .true.
    return
  else
    Is_First_Iter = .false.
    return
  end if
end if

! If Slapaf information has been stripped out (e.g. IRC restart)
call qpg_iArray('Slapaf Info 1',Found,Tmp(1))
if (Found) then
  call Get_iArray('Slapaf Info 1',Tmp,7)
  if (Tmp(1) == -99) then
    Is_First_Iter = .true.
    return
  end if
end if

! If MOLCAS_ITER <= 1
call Getenvf('MOLCAS_ITER',EnvStr)
read(EnvStr,*) Iter
if (Iter <= 1) then
  Is_First_Iter = .true.
  return
end if

Is_First_Iter = .false.

end function Is_First_Iter
