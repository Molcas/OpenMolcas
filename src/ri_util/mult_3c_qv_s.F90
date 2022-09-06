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

subroutine Mult_3C_Qv_s(A_3C,nA_3C,Qv,nQv,Rv,n_Rv,nVec,iOff_3C,nIrrep,Out_of_Core,Lu_Q,QMode)
!***********************************************************************
!     Author:   F. Aquilante                                           *
!                                                                      *
!     Qv: is a symmetry blocked square matrix                          *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
real*8 A_3C(nA_3C), Qv(nQv), Rv(n_Rv)
integer iOff_3C(3,0:nIrrep-1), nVec(0:7), Lu_Q(0:nIrrep-1)
logical Out_of_Core
character*1 QMode

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

if (Out_of_Core) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iIrrep=0,nIrrep-1
    nI = iOff_3C(2,iIrrep)
    nMuNu = iOff_3C(1,iIrrep)
    if ((nMuNu <= 0) .or. (nI <= 0)) Go To 999

    iOffR2 = iOffR
    iOffA2 = iOffA
    mQv = nI*nVec(iIrrep)
    iAddr = 0
    do while (mQv >= nI)

      nK = min(mQv,nQv)/nI
      lQv = nI*nK
      call dDaFile(Lu_Q(iIrrep),2,Qv,lQv,iAddr)

      call A_3C_Qv_s(A_3C(iOffA2),Qv,Rv(iOffR2),nMuNu,nI,nK,QMode)
      mQv = mQv-lQv
      iOffR2 = iOffR2+lstepR*nMuNu*nK
      iOffA2 = iOffA2+lstepA*nMuNu*nK
    end do

    iOffA = iOffA+nMuNu*nI
    iOffR = iOffR+nMuNu*nVec(iIrrep)
999 continue
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else  ! In-Core
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  iOffQ = 1
  do iIrrep=0,nIrrep-1
    nI = iOff_3C(2,iIrrep)
    nMuNu = iOff_3C(1,iIrrep)
    if ((nMuNu <= 0) .or. (nI <= 0)) Go To 998

    call A_3C_Qv_s(A_3C(iOffA),Qv(iOffQ),Rv(iOffR),nMuNu,nI,nVec(iIrrep),QMode)
998 continue
    iOffA = iOffA+nMuNu*nI
    iOffR = iOffR+nMuNu*nVec(iIrrep)
    iOffQ = iOffQ+nI*nVec(iIrrep)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mult_3C_Qv_s
