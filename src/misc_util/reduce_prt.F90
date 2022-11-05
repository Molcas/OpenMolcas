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

function Reduce_Prt()

use UnixInfo, only: ProgName, SuperName
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: Reduce_Prt
integer(kind=iwp) :: i, Err
character(len=80) :: Word

Reduce_Prt = .false.

! Do not reduce printing in last_energy

if (SuperName == 'last_energy') return

! Reduce printing if iter > 1

call GetEnvF('MOLCAS_ITER',Word)
read(Word,*) i
if (i > 1) Reduce_Prt = .true.

! ... but not if MOLCAS_REDUCE_PRT = NO

if (Reduce_Prt) then
  call GetEnvF('MOLCAS_REDUCE_PRT',Word)
  if (Word(1:1) == 'N') Reduce_Prt = .false.
end if

! ... or if we are not inside a loop (EMIL_InLoop < 1)

if (Reduce_Prt) then
  call GetEnvF('EMIL_InLoop',Word)
  i = 0
  read(Word,*,IOSTAT=Err) i
  if (i < 1) Reduce_Prt = .false.
end if

! ... or if first iteration of a saddle branch (SADDLE_FIRST = 1)

if (Reduce_Prt) then
  call GetEnvF('SADDLE_FIRST',Word)
  i = 0
  read(Word,*,IOSTAT=Err) i
  if (i == 1) Reduce_Prt = .false.
end if

! In any case, reduce printing inside numerical gradients,
! unless specified otherwise (MOLCAS_REDUCE_NG_PRT = NO).

if (.not. Reduce_Prt) then
  if ((SuperName == 'numerical_gradient') .and. (ProgName /= 'numerical_gradient')) then
    call GetEnvF('MOLCAS_REDUCE_NG_PRT',Word)
    if (Word(1:1) /= 'N') Reduce_Prt = .true.
  end if
end if

return

end function Reduce_Prt
