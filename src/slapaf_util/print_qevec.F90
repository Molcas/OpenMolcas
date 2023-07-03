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

subroutine Print_qEVec(EVec,nH,EVal,nq,rK,qEVec,LuTmp)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 EVec(nH,nH), rK(nq,nH), qEVec(nq,nH), EVal(nH*(nH+1)/2)
character(len=14) qLbl(nq)

do iq=1,nq
  read(LuTmp) qLbl(iq),(rK(iq,iQQ),iQQ=1,nH)
end do

call DGEMM_('N','N',nq,nH,nH,1.0d0,rK,nq,EVec,nH,0.0d0,qEVec,nq)

Lu = 6
Thr = 0.0001d0
IncQQ = 5
do iiQQ=1,nH,IncQQ
  mQQ = min(nH,iiQQ+IncQQ-1)
  write(Lu,*)
  write(Lu,'(14X,5I10)') (iQQ,iQQ=iiQQ,mQQ)
  write(Lu,'(1X,A,5F10.6)') 'Eigenvalues   ',(EVal(iQQ*(iQQ+1)/2),iQQ=iiQQ,mQQ)
  write(Lu,*)
  do iq=1,nq
    temp = sqrt(DDot_(nH,qEVec(iq,1),nq,qEVec(iq,1),nq))
    if (temp > Thr) write(Lu,'(1X,A,5F10.6)') qLbl(iq),(qEVec(iq,iQQ),iQQ=iiQQ,mQQ)
  end do
  write(Lu,*)
end do

return

end subroutine Print_qEVec
