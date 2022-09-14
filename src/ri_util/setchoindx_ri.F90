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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine SetChoIndx_RI(iiBstRSh,nnBstRSh,IndRed,IndRsh,iRS2F,I_nSym,I_nnShl,I_mmBstRT,iShij,nShij)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use ChoArr, only: iSP2F, iBasSh, nBasSh, nBstSh
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: I_nSym, I_nnShl, I_mmBstRT, iRS2F(2,I_mmBstRT), nShij, iShij(2,nShij)
integer(kind=iwp), intent(out) :: iiBstRSh(I_nSym,I_nnShl,3), nnBstRSh(I_nSym,I_nnShl,3), IndRed(I_mmBstRT,3), IndRsh(I_mmBstRT)
#include "cholesky.fh"
integer(kind=iwp) :: i, ia, iaa, iab, ib, ibb, iCount, iRS(8), iSh_ij, iShla, iShlab, iShlb, iSym, iSyma, iSymb, jRS, jRS1, jRS2, &
                     nErr
integer(kind=iwp), external :: Cho_iSAOSh

! nnBstRSh(iSym,iSh_ij,1) = #elements in compound sym. iSym of
!                           shell-pair ab in 1st reduced set.
! IndRSh(jRS): shell-pair to which element jRS of first reduced set
!              belongs.
! IndRed(jRS,1): address (without symmetry) in shell-pair of element
!                jRS of first reduced set.
! ------------------------------------------------------------------

nnBstRSh(:,:,1) = 0
iRS(1:nSym) = iiBstR(1:nSym,1)
do iSh_ij=1,nShij
  iShla = iShij(1,iSh_ij)
  iShlb = iShij(2,iSh_ij)
  iShlab = iTri(iShla,iShlb)
  !write(u6,*) 'iSh_ij,iShlab,iShla,iShlb=',iSh_ij,iShlab,iShla,iShlb
  if (iShlab /= iSP2F(iSh_ij)) call SysAbendMsg('SetChoIndx_RI','SP2F setup error',' ')

  if (iShla > iShlb) then

    ! code for shell a > shell b

    do iSymb=1,nSym
      do ibb=1,nBasSh(iSymb,iShlb)
        ib = iBasSh(iSymb,iShlb)+ibb
        do iSyma=1,nSym
          iSym = Mul(iSyma,iSymb)
          do iaa=1,nBasSh(iSyma,iShla)
            ia = iBasSh(iSyma,iShla)+iaa
            iab = nBstSh(iShla)*(ib-1)+ia
            nnBstRSh(iSym,iSh_ij,1) = nnBstRSh(iSym,iSh_ij,1)+1
            iRS(iSym) = iRS(iSym)+1
            IndRSh(iRS(iSym)) = iShlab
            IndRed(iRS(iSym),1) = iab
          end do
        end do
      end do
    end do

  else

    ! code for shell a = shell b follows

    do ia=1,nBstSh(iShla)
      iSyma = Cho_iSAOSh(ia,iShla)
      do ib=1,ia
        iab = iTri(ia,ib)
        iSymb = Cho_iSAOSh(ib,iShlb)
        iSym = Mul(iSyma,iSymb)
        nnBstRSh(iSym,iSh_ij,1) = nnBstRSh(iSym,iSh_ij,1)+1
        iRS(iSym) = iRS(iSym)+1
        IndRSh(iRS(iSym)) = iShlab
        IndRed(iRS(iSym),1) = iab
      end do
    end do

  end if
end do   ! iSh_ij

! Check.
! ------

nErr = 0
do iSym=1,nSym
  iCount = nnBstRSh(iSym,1,1)
  do iSh_ij=2,nnShl
    iCount = iCount+nnBstRSh(iSym,iSh_ij,1)
  end do
  if (iCount /= nnBstR(iSym,1)) then
    nErr = nErr+1
  end if
end do
if (nErr /= 0) then
  call SysAbendMsg('SetChoIndx_RI','Setup error','iCount vs. nnBstR')
end if
do iSym=1,nSym
  if ((iRS(iSym)-iiBstR(iSym,1)) /= nnBstR(iSym,1)) nErr = nErr+1
end do
if (nErr /= 0) then
  call SysAbendMsg('SetChoIndx_RI','Setup error','ShP RS1 count')
end if

! iiBstRSh(iSym,iSh_ij,1) = offset to elements in compound sym. iSym
!                           of shell-pair ab in 1st reduced set.
! ------------------------------------------------------------------

do iSym=1,nSym
  iiBstRSh(iSym,1,1) = 0
  do iSh_ij=2,nnShl
    iiBstRSh(iSym,iSh_ij,1) = iiBstRSh(iSym,iSh_ij-1,1)+nnBstRSh(iSym,iSh_ij-1,1)
  end do
end do

! Check.
! ------

nErr = 0
do iSym=1,nSym
  do iSh_ij=1,nnShl
    jRS1 = iiBstR(iSym,1)+iiBstRSh(iSym,iSh_ij,1)+1
    jRS2 = jRS1+nnBstRSh(iSym,iSh_ij,1)-1
    do jRS=jRS1,jRS2
      if (IndRSh(jRS) /= iSP2F(iSh_ij)) nErr = nErr+1
    end do
  end do
end do
if (nErr /= 0) then
  call SysAbendMsg('SetChoIndx_RI','Setup error','IndRSh')
end if

! Copy index arrays to "locations" 2 and 3.
! Note: IndRed here returns the index in 1st reduced set.
! -------------------------------------------------------

do i=2,3
  do jRS=1,nnBstRT(1)
    IndRed(jRS,i) = jRS
  end do
  iiBstRSh(:,:,i) = iiBstRSh(:,:,1)
  nnBstRSh(:,:,i) = nnBstRSh(:,:,1)
end do

call Cho_RStoF(iRS2F,2,nnBstRT(1),1)

return

end subroutine SetChoIndx_RI
