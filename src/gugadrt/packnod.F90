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

implicit real*8(a-h,o-z)
dimension ibuf(lbuf)
integer*4, parameter :: one4 = 1
integer, parameter :: i4 = kind(one4)

inv = ival
nimod = mod(idx,nin)
if (nimod == 0) then
  ngrp = idx/nin
  nidbit = 0
else
  ngrp = idx/nin+1
  nidbit = (nin-nimod)*nbit
end if

!write(6,"(b64.64)") inv
#ifdef _AIX_
call abend()
#else
!call mvbits(inv,0,nbit,ibuf(ngrp),nidbit)
call mvbits(inv,0,int(nbit,i4),ibuf(ngrp),int(nidbit,i4))
#endif
!write(6,"(b64.64)") ibuf(ngrp)

return

end subroutine packnod
