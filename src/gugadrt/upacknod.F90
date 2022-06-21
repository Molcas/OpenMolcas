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

subroutine upacknod(ibuf,idx,ival,nin,nbit,lbuf)

!**************************************************
! pack integral ival into ibuf on bit mode
! ibuf() integral buffer array
! idx    index
! nin    number of integrals in one integral
! lbuf   length of ibuf

use Definitions, only: iwp, DefInt

implicit none
integer(kind=iwp), intent(in) :: lbuf, ibuf(lbuf), idx, nin, nbit
integer(kind=iwp), intent(out) :: ival
integer(kind=iwp) :: isp, ngrp, nidbit, nimod

nimod = mod(idx,nin)
if (nimod == 0) then
  ngrp = idx/nin
  nidbit = 0
else
  ngrp = idx/nin+1
  nidbit = (nin-nimod)*nbit
end if

isp = nidbit
!inv = ibuf(ngrp)
!ival = ibits(inv,isp,nbit-1)
ival = 0
#ifdef _AIX_
call abend()
#else

!call mvbits(ibuf(ngrp),isp,nbit,ival,0)
call mvbits(ibuf(ngrp),int(isp,kind=DefInt),int(nbit,kind=DefInt),ival,0)
#endif
!write(u6,*) isp,nbit-1
!write(u6,'(b64.64)') ival

return

end subroutine upacknod
