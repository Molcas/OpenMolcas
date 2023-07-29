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
!               2021, Roland Lindh                                     *
!***********************************************************************

subroutine CHO_GetShFull(LabJ,lLabJ,JNUM,JSYM,IREDC,ChoV,SvShp,mmShl,iShp_rs,mmShl_tot)

use Symmetry_Info, only: Mul
use Cholesky, only: iBasSh, iiBstR, IndRed, IndRSh, iRS2F, iShlSO, iSOShl, nDimRS, nnBstR, nnShl_tot
use Cholesky_Structures, only: L_Full_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: lLabJ, JNUM, JSYM, IREDC, mmShl, mmShl_tot, iShp_rs(mmShl_tot)
real(kind=wp), intent(in) :: LabJ(lLabJ)
type(L_Full_Type), intent(_OUT_) :: ChoV
real(kind=wp), intent(out) :: SvShp(mmShl,2)
integer(kind=iwp) :: i1, iag, ias, iaSg, iaSh, ibg, ibs, ibSg, ibSh, iLoc, iRab, iShp, iSyma, iSymb, jRab, jShp, JVEC, kLabJ, &
                     kRab, NREAD
integer(kind=iwp), external :: cho_isao

!*********************************************************
!
! From Reduced sets to full storage
! ---------------------------------
!
!  L{a,b,J} ---> L(a,J,b)
!
!*********************************************************

iLoc = 3 ! use scratch location in reduced index arrays

ChoV%A0(:) = Zero
SvShp(:,:) = Zero

if (JSYM == 1) then

  NREAD = 0

  do JVEC=1,JNUM

    kLabJ = NREAD
    NREAD = NREAD+nDimRS(jSym,IREDC)

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,iLoc) ! addr in 1st red set

      iShp = IndRSh(iRab) ! shell pair to which it belongs

      iag = iRS2F(1,iRab) ! global address
      ibg = iRS2F(2,iRab)

      iaSh = iSOShl(iag) ! shell to which it belongs
      ibSh = iSOShl(ibg) ! iaSh >= ibSh !!!!!!

      iaSg = iShlSO(iag) ! index of SO within its shell
      ibSg = iShlSO(ibg)

      iSyma = cho_isao(iag) ! symmetry block sym(a)=sym(b)

      ias = iaSg-iBasSh(iSyma,iaSh) ! addr within its shell
      ibs = ibSg-iBasSh(isyma,ibSh)

      kLabJ = kLabJ+1

      ChoV%SPB(iSyma,iShp_rs(iShp),1)%A3(ias,JVEC,ibs) = LabJ(kLabJ)

      i1 = 1
      if (ibSh /= iaSh) i1 = 2

      ChoV%SPB(iSyma,iShp_rs(iShp),i1)%A3(ibs,JVEC,ias) = LabJ(kLabJ)

      SvShp(iShp_rs(iShp),2) = SvShp(iShp_rs(iShp),2)+LabJ(kLabJ)**2

    end do

    do jShp=1,nnShl_tot ! Maximize over vectors
      if (iShp_rs(jShp) > 0) then
        SvShp(iShp_rs(jShp),1) = max(SvShp(iShp_rs(jShp),1),SvShp(iShp_rs(jShp),2))
        SvShp(iShp_rs(jShp),2) = zero
      end if
    end do

  end do

else

  NREAD = 0

  do JVEC=1,JNUM

    kLabJ = NREAD
    NREAD = NREAD+nDimRS(jSym,IREDC)

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,iLoc) ! addr in 1st red set

      iShp = IndRSh(iRab) ! shell pair to which it belongs

      iag = iRS2F(1,iRab) ! global address
      ibg = iRS2F(2,iRab)

      iaSh = iSOShl(iag) ! shell to which it belongs
      ibSh = iSOShl(ibg) ! ibsh<=>iaSh

      iaSg = iShlSO(iag) ! index of SO within its shell
      ibSg = iShlSO(ibg)

      iSyma = cho_isao(iag) ! symmetry block
      iSymb = Mul(jSym,iSyma) ! iSyma >= iSymb

      i1 = 1
      if (iaSh < ibSh) i1 = 2

      ias = iaSg-iBasSh(iSyma,iaSh) ! addr within its shell
      ibs = ibSg-iBasSh(iSymb,ibSh)

      kLabJ = kLabJ+1

      ChoV%SPB(iSyma,iShp_rs(iShp),i1)%A3(ias,JVEC,ibs) = LabJ(kLabJ)

      SvShp(iShp_rs(iShp),2) = SvShp(iShp_rs(iShp),2)+LabJ(kLabJ)**2

    end do

    do jShp=1,nnShl_tot
      if (iShp_rs(jShp) > 0) then
        SvShp(iShp_rs(jShp),1) = max(SvShp(iShp_rs(jShp),1),SvShp(iShp_rs(jShp),2))
        SvShp(iShp_rs(jShp),2) = zero
      end if
    end do

  end do

end if

return

end subroutine CHO_GetShFull
