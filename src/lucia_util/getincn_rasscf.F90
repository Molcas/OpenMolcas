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

subroutine GETINCN_RASSCF(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,IPNT2,NSMOB,INH1,ICOUL)
! interface to RASSCF common blocks

use wadr, only: tuvx

implicit none
integer ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, IKSM, JLSM, NSMOB, ICOUL
! For Jesper and openMP
integer IPNT2(*), INH1(*)
real*8 XINT(*)

call GETINCN_RASSCFS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,TUVX,NSMOB,ICOUL)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(IPNT2)
  call Unused_integer_array(INH1)
end if

end subroutine GETINCN_RASSCF
