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

subroutine Mult_Vk_Qv_s(V_k,nV_k,Qv,nQv,V_kQ,nV_kQ,nBas_Aux,nVec,nIrrep,QMode)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nV_k, nQv, nV_kQ, nIrrep, nBas_Aux(0:nIrrep-1), nVec(0:nIrrep-1)
real(kind=wp), intent(in) :: V_k(nV_k)
real(kind=wp), intent(out) :: Qv(nQv), V_kQ(nV_kQ)
character, intent(in) :: QMode
integer(kind=iwp) :: iAddr, iIrrep, iOffA, iOffB, iSeed, kp_V_k, lQv, lstepA, lstepB, Lu_Q, mQv, nI, nJ, nK, nMuNu
logical(kind=iwp) :: Out_of_Core
character(len=6) :: Name_Q
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
nMuNu = 1
kp_V_k = 1
lstepA = 0
lstepB = 1
if (Qmode == 'T') then
  V_kQ(:) = Zero
  lstepA = 1
  lstepB = 0
end if

do iIrrep=0,0  ! loop is wisely restricted to tot. symm. irrep
  iSeed = 55+iIrrep
  Lu_Q = IsFreeUnit(iSeed)
  write(Name_Q,'(A4,I2.2)') 'QVEC',iIrrep
  call DaName_MF_WA(Lu_Q,Name_Q)
  iAddr = 0

  nI = nBas_Aux(iIrrep)
  if (iIrrep == 0) nI = nI-1
  nJ = nVec(iIrrep)
  mQv = nI*nJ
  Out_of_Core = mQv > nQv
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (.not. Out_of_Core) then

    ! in-core case

    call dDaFile(Lu_Q,2,Qv,mQv,iAddr)

    call A_3C_Qv_s(V_k(kp_V_k),Qv,V_kQ,nMuNu,nI,nJ,QMode)

  else

    ! Out-of-core case

    iOffA = kp_V_k
    iOffB = iOffA
    do while (mQv >= nI)

      nK = min(mQv,nQv)/nI
      lQv = nI*nK
      call dDaFile(Lu_Q,2,Qv,lQv,iAddr)

      call A_3C_Qv_s(V_k(iOffA),Qv,V_kQ(iOffB),nMuNu,nI,nK,QMode)
      mQv = mQv-lQv
      iOffA = iOffA+lstepA*nK
      iOffB = iOffB+lstepB*nK
    end do

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call DaClos(Lu_Q)

  kp_V_k = kp_V_k+nBas_Aux(iIrrep)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mult_Vk_Qv_s
