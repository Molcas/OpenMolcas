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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************

subroutine Find_Basis_Set(DirName,ExtBasDir,bType)
!***********************************************************************
! Object: make a guess for a file                                      *
! Author: Valera Veryazov, Theoretical Chemistry, Chemical Centre      *
!         University of Lund, Lund, Sweden                             *
!***********************************************************************

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: DirName
character(len=*), intent(in) :: ExtBasDir, bType
integer(kind=iwp) :: i, iAbsName
logical(kind=iwp) :: Exists
character(len=512) :: tmp
character(len=256) :: Molcas, CurrDir

if (ExtBasDir /= ' ') then
  i = index(ExtBasDir,' ')
  iAbsName = 0
  if (ExtBasDir(1:1) /= '/') then
    iAbsName = 1
    CurrDir = ' '
    call getenvf('CurrDir',CurrDir)
  end if
  if (iAbsName == 0) then
    tmp = ExtBasDir(1:i-1)//'/'//bType
  else
    tmp = CurrDir(1:index(CurrDir,' ')-1)//'/'//ExtBasDir(1:i-1)//'/'//bType
  end if
  !write(u6,*) tmp
  call f_Inquire(tmp,Exists)
  if (Exists) then
    if (iAbsName == 0) then
      tmp = ExtBasDir(1:i-1)
    else
      tmp = CurrDir(1:index(CurrDir,' ')-1)//'/'//ExtBasDir(1:i-1)
    end if
    i = index(tmp,' ')
    DirName = tmp(1:i-1)
    !write(u6,*) '>',DirName(1:i-1),'<'
    return
  end if
end if
if (DirName(1:13) /= 'basis_library') return

Molcas = ' '
call getenvf('MOLCAS_BASIS',Molcas)
if (Molcas == ' ') then
  call getenvf('MOLCAS',Molcas)
  i = index(Molcas,' ')
  DirName = Molcas(1:i-1)//'/basis_library'
else
  i = index(Molcas,' ')
  DirName = Molcas(1:i-1)
end if
i = index(DirName,' ')
if (i == 0) then
  call WarningMessage(2,'Too long path to Molcas')
  call Abend()
end if

return

end subroutine Find_Basis_Set
