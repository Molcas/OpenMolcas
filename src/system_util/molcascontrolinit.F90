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
! Copyright (C) 2000-2016, Valera Veryazov                             *
!***********************************************************************

subroutine MolcasControlInit(label)

parameter(nLines=20)
character*(*) Label
character*512 tmp
character*32 My
character filename*16
integer StrnLn

filename = 'molcas.control'
iC = 0
tmp = Label(1:len(Label))
Lu = 1
Lu = isfreeunit(Lu)
open(Lu,File=filename)
write(Lu,'(a)') '# Molcas control file: change # to ! to activate.'
islast = 0

10 i = index(tmp,',')
My = ' '
if (i > 0) then
  My(1:i-1) = tmp(1:i-1)
  tmp = tmp(i+1:)
else
  My = trim(tmp)
  islast = 1
end if
iC = iC+1
if (ic > nLines) call abend()
i = StrnLn(My)
if (index(My,'=') == 0) then
  i = i+1
  My(i:i) = '='
end if
write(Lu,'(a,a)') '#',My(1:i)
if (islast == 0) goto 10
close(Lu)

return

end subroutine MolcasControlInit
