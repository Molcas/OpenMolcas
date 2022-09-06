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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine Mult_Zp_Qv_s(Zp,nZp,Qv,nQv,Rv,n_Rv,nVec,nMuNu,nI,nIrrep,QMode)
!***********************************************************************
!     Author:   F. Aquilante                                           *
!                                                                      *
!     Qv: is a symmetry blocked square matrix                          *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
real*8 Zp(nZp), Qv(nQv), Rv(n_Rv)
integer nVec(0:nIrrep-1), nMuNu(0:nIrrep-1), nI(0:nIrrep-1)
character QMode*1, Name_Q*6

!                                                                      *
!***********************************************************************
!                                                                      *
lstepA = 0
lstepR = 1
if (Qmode == 'T') then
  call FZero(Rv,n_Rv)
  lstepA = 1
  lstepR = 0
end if

iOffA = 1
iOffR = 1

do iIrrep=0,nIrrep-1
  nI_ = nI(iIrrep)
  if (iIrrep == 0) nI_ = nI_-1
  nMuNu_ = nMuNu(iIrrep)
  if ((nMuNu_ <= 0) .or. (nI_ <= 0)) Go To 999

  iSeed = 55+iIrrep
  Lu_Q = IsFreeUnit(iSeed)
  write(Name_Q,'(A4,I2.2)') 'QVEC',iIrrep
  call DaName_MF_WA(Lu_Q,Name_Q)
  iAddr = 0

  iOffR2 = iOffR
  iOffA2 = iOffA
  mQv = nI_*nVec(iIrrep)
  do while (mQv >= nI_)

    nK = min(mQv,nQv)/nI_
    lQv = nI_*nK
    call dDaFile(Lu_Q,2,Qv,lQv,iAddr)

    call A_3C_Qv_s(Zp(iOffA2),Qv,Rv(iOffR2),nMuNu_,nI_,nK,QMode)
    mQv = mQv-lQv
    iOffR2 = iOffR2+lstepR*nMuNu_*nK
    iOffA2 = iOffA2+lstepA*nMuNu_*nK
  end do

  iOffA = iOffA+nMuNu_*nI_
  iOffR = iOffR+nMuNu_*nVec(iIrrep)
  call DaClos(Lu_Q)
999 continue
end do

return

end subroutine Mult_Zp_Qv_s
