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

subroutine packnod(ibuf,idx,ival,nin,nbit,lbuf)
!**************************************************
! pack integral ival into ibuf on bit mode
! ibuf() integral buffer array
! idx    index
! nin    number of integrals in one integral
! lbuf   length of ibuf

use Definitions, only: iwp, DefInt

implicit none
integer(kind=iwp), intent(in) :: idx, ival, nin, nbit, lbuf
integer(kind=iwp), intent(inout) :: ibuf(lbuf)
integer(kind=iwp) :: inv, ngrp, nidbit, nimod

inv = ival
nimod = mod(idx,nin)
if (nimod == 0) then
  ngrp = idx/nin
  nidbit = 0
else
  ngrp = idx/nin+1
  nidbit = (nin-nimod)*nbit
end if

!write(u6,'(b64.64)') inv
#ifdef _AIX_
call abend()
#else
!call mvbits(inv,0,nbit,ibuf(ngrp),nidbit)
call mvbits(inv,0,int(nbit,kind=DefInt),ibuf(ngrp),int(nidbit,kind=DefInt))
#endif
!write(u6,'(b64.64)') ibuf(ngrp)

return

end subroutine packnod
