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
!***********************************************************************
!                                                                      *
! Author:   Valera Veryazov 2000-2016                                  *
!           Theoretical Chemistry                                      *
!           Lund University                                            *
!           Sweden                                                     *
!                                                                      *
!***********************************************************************
!  MolcasControl
!
!> @brief
!>   Query a string from the control file
!> @author V. Veryazov
!>
!> @details
!> Only lines started from ``!`` are in use.
!>
!> If user modified molcas.control file
!> (by placing ``!`` instead of ``#``)
!> return a value of the field (as a string)
!> and mark the label as a comment,
!> else return blank value.
!>
!> Usage:
!>
!> \code
!> Call MolcasControl('SHUTDOWN',Value)
!> if(Value.eq.'YES') Call abend
!> \endcode
!>
!> @side_effects
!> file molcas.control
!>
!> @param[in]  Label Query string
!> @param[out] Value Returned value
!***********************************************************************

subroutine MolcasControl(Label,value)

parameter(nLines=20)
character filename*16
character*80 Line(nLines)
character*(*) Label
character*(*) value
logical Exist
logical Modify
integer StrnLn

filename = 'molcas.control'
value = ' '
Modify = .false.
value = ' '
call f_inquire(filename,Exist)
if (.not. Exist) return
Lu = 1
Lu = isfreeunit(Lu)
open(Lu,File=filename)
iLine = 1
1 read(Lu,'(a)',end=100,err=100) Line(iLine)
if (Line(iLine)(1:1) == '!') Modify = .true.
iLine = iLine+1
if (iLine < nLines) goto 1

100 continue
close(Lu)

if (.not. Modify) return

open(Lu,File=filename)
do ic=1,iLine-1
  if (Line(ic)(1:1) == '!') then
    i = index(Line(ic)(2:),'=')
    if (i > 0) then
      if (Line(ic)(2:i) == Label) then
        Line(ic)(1:1) = '#'
        Modify = .true.
        value = Line(ic)(i+2:)
      end if
    end if
  end if
  i = StrnLn(Line(ic))
  write(Lu,'(a)') Line(ic)(1:i)
end do
close(Lu)

return

end subroutine MolcasControl
