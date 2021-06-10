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

subroutine prgmtranslate(in,out,lout)

character*(*) in, out
#ifdef _DEBUGPRINT_IO_
character LL*128
#endif
integer Strnln
external Strnln

lin = Strnln(in)
out = ' '
if (index(in,'/') /= 0) then
  ! just in case if we processing translated name!
  out = in
  lout = lin
else
  call prgmtranslatec(in,lin,out,lout,1)
end if
out = out(1:lout)
#ifdef _DEBUGPRINT_IO_
LL = 'QWERTYUIOPASDFGHJKLZXCVBNM1234567890qwertyuiopasdfghjklzxcvbnm /.-_*'

!write(6,*) 'Translate: >',in(1:lin),'< to >',out(1:lout),'<'
do i=1,lin
  if (index(LL,in(i:i)) == 0) then
    write(6,*) 'Translate: >',in(1:lin),'< to >',out(1:lout),'<'
    write(6,*) 'invalid in:',in(i:i),'<'
    call abend()
  end if
end do
do i=1,lout
  if (index(LL,out(i:i)) == 0) then
    write(6,*) 'Translate: >',in(1:lin),'< to >',out(1:lout),'<'
    write(6,*) 'invalid out:',out(i:i),'<'
    call abend()
  end if
end do
#endif

return

end subroutine prgmtranslate
