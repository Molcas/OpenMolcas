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
! Copyright (C) 2002, Roland Lindh                                     *
!***********************************************************************

subroutine Cmbn_EF_DPnt(EF,nTs,DPnt,MxAto,DCntr,nS,iSph,Q,Grad,nGrad)
!***********************************************************************
!                                                                      *
!  Combine EF with DPnt array.                                         *
!                                                                      *
!  Roland Lindh                                                        *
!  020117                                                              *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use Disp, only: IndDsp
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTs, MxAto, nS, iSph(nTs), nGrad
real(kind=wp), intent(in) :: EF(3,2,nTs), DPnt(nTs,MxAto,3,3), DCntr(nS,MxAto,3,3), Q(2,nTs)
real(kind=wp), intent(inout) :: Grad(nGrad)
integer(kind=iwp) :: iCar, iCen, iCnt, iCnttp, iComp, iIrrep, iTs, jSph, mdc, nDispS
real(kind=wp) :: QTot
real(kind=wp), parameter :: tol = 1.0e-8_wp
logical(kind=iwp), external :: TF

iIrrep = 0

mdc = 0
iCen = 1
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux) cycle
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    nDispS = IndDsp(mdc,iIrrep)

    do iCar=0,2
      iComp = 2**iCar
      if (TF(mdc,iIrrep,iComp)) then
        nDispS = nDispS+1
        do iTs=1,nTs
          jSph = iSph(iTs)
          QTot = Q(1,iTs)+Q(2,iTs)
          Grad(nDispS) = Grad(nDispS)+QTot*((EF(1,1,iTs)-EF(1,2,iTs))*(DPnt(iTs,iCen,iCar+1,1)+DCntr(jSph,iCen,iCar+1,1))+ &
                                            (EF(2,1,iTs)-EF(2,2,iTs))*(DPnt(iTs,iCen,iCar+1,2)+DCntr(jSph,iCen,iCar+1,2))+ &
                                            (EF(3,1,iTs)-EF(3,2,iTs))*(DPnt(iTs,iCen,iCar+1,3)+DCntr(jSph,iCen,iCar+1,3)))
        end do
      end if
    end do
    iCen = iCen+nIrrep/dc(mdc)%nStab
  end do
end do

return

end subroutine Cmbn_EF_DPnt
