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

subroutine Integral_WrOut_Cho_diag( &
#                                  define _CALLING_
#                                  include "int_wrout_interface.fh"
                                  )
! calls the proper routines IndSft/PLF
!    if IntOrd_jikl==.TRUE. integral order within symblk: jikl
!                     else  integral order within symblk: ijkl

implicit real*8(A-H,O-Z)
# define _FIXED_FORMAT_
#include "int_wrout_interface.fh"

! call sorting routine

if (mSym == 1) then
  call PLF_Cho_Diag(TInt,nTInt,AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,Shijij .and. IJeqKL,iBas,jBas,kBas, &
                    lBas,kOp)
else
  call IndSft_Cho_Diag(TInt,nTInt,iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
end if

return
if (.false.) then
  call Unused_integer(nSkal)
end if

end subroutine Integral_WrOut_Cho_diag
