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
! Copyright (C) 2001-2016, Valera Veryazov                             *
!***********************************************************************

subroutine prgmtranslate(namein,nameout,lout)

use Definitions, only: iwp
#ifdef _DEBUGPRINT_IO_
use Definitions, only: u6
#endif

implicit none
character(len=*), intent(in) :: namein
character(len=*), intent(out) :: nameout
integer(kind=iwp), intent(out) :: lout
integer(kind=iwp) :: lin
#ifdef _DEBUGPRINT_IO_
character(len=68), parameter :: LL = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz /.-_*'
#endif
integer(kind=iwp), external :: Strnln

lin = Strnln(namein)
nameout = ' '
if (index(namein,'/') /= 0) then
  ! just in case we process a translated name!
  nameout = namein
  lout = lin
else
  call prgmtranslatec(namein,lin,nameout,lout,1)
end if
nameout = nameout(1:lout)
#ifdef _DEBUGPRINT_IO_
!write(u6,*) 'Translate: >',namein(1:lin),'< to >',nameout(1:lout),'<'
do i=1,lin
  if (index(LL,namein(i:i)) == 0) then
    write(u6,*) 'Translate: >',namein(1:lin),'< to >',nameout(1:lout),'<'
    write(u6,*) 'invalid namein:',namein(i:i),'<'
    call abend()
  end if
end do
do i=1,lout
  if (index(LL,nameout(i:i)) == 0) then
    write(u6,*) 'Translate: >',namein(1:lin),'< to >',nameout(1:lout),'<'
    write(u6,*) 'invalid nameout:',nameout(i:i),'<'
    call abend()
  end if
end do
#endif

return

end subroutine prgmtranslate
