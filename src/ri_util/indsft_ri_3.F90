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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine IndSft_RI_3(iCmp,iShell,jBas,kBas,lBas,iAO,iAOst,ijkl,SOint,nSOint,TInt,nTInt,iOff,iShlSO,nBasSh,iSOShl,nSO, &
                       nShell,nSym,iSSOff)
!***********************************************************************
!  object: to sift and index the SO integrals.                         *
!                                                                      *
!          the indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          april '90                                                   *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use Symmetry_Info, only: Mul, nIrrep
use sort_data, only: nSkip
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iShell(4), jBas, kBas, lBas, iAO(4), iAOst(4), ijkl, nSOint, nTInt, iOff(3,0:7), nSO, &
                                 iShlSO(nSO), nShell, nSym, nBasSh(0:nSym-1,nShell), iSOShl(nSO), iSSOff(0:nIrrep-1,0:nIrrep-1)
real(kind=wp), intent(in) :: SOint(ijkl,nSOint)
real(kind=wp), intent(inout) :: TInt(nTInt)
integer(kind=iwp) :: i2, i3, i4, iAux, iC, iD, iIrrep, iOff_L, iShC, iShD, ix, j, j1, j12, j2, j3, j4, jOffSO(0:7), jSO, jSOj, &
                     jSym(0:7), kl, kl_B, kSO, kSOk, kSym(0:7), lCmpMx, lSO, lSOl, lSym(0:7), memSO2, mm, n3C, nC, nD, nijkl
logical(kind=iwp) :: qkl, Shkl

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
jOffSO(0) = 0
do iIrrep=1,nIrrep-1
  jOffSO(iIrrep) = jOffSO(iIrrep-1)+nBas(iIrrep-1)
end do
memSO2 = 0

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.

Shkl = iShell(3) == iShell(4)
if (iShell(4) > iShell(3)) then
  call WarningMessage(2,'Error in IndSft_RI_3')
  write(u6,*) 'iShell(4) > iShell(3)'
  call Abend()
end if

j1 = 0
do i2=1,iCmp(2)
  do j=0,nIrrep-1
    ix = 0
    if (iAOtSO(iAO(2)+i2,j) > 0) ix = 2**j
    jSym(j) = ix
  end do
  do i3=1,iCmp(3)
    do j=0,nIrrep-1
      ix = 0
      if (iAOtSO(iAO(3)+i3,j) > 0) ix = 2**j
      kSym(j) = ix
    end do
    lCmpMx = iCmp(4)
    if (Shkl) lCmpMx = i3
    do i4=1,lCmpMx
      do j=0,nIrrep-1
        ix = 0
        if (iAOtSO(iAO(4)+i4,j) > 0) ix = 2**j
        lSym(j) = ix
      end do
      qkl = i3 == i4

      ! loop over Irreps which are spanned by the basis function.
      ! again, the loop structure is restricted to ensure unique
      ! integrals.

      do j2=0,nIrrep-1
        if (jSym(j2) == 0) cycle
        j12 = Mul(j1+1,j2+1)-1

        do j3=0,nIrrep-1
          if (kSym(j3) == 0) cycle
          j4 = Mul(j12+1,j3+1)-1
          if (lSym(j4) == 0) cycle
          if (Shkl .and. qkl .and. (j4 > j3)) cycle

          memSO2 = memSO2+1
          if ((nSkip(j2+1)+nSkip(j3+1)+nSkip(j4+1)) /= 0) cycle
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Number of auxiliary basis functions in this symmetry block.
          mm = iOff(1,j12)
          if (mm == 0) cycle
          ! Effective number of valence basis products in this symmetry block.
          n3C = iOff(3,j12)
          if (n3C == 0) cycle
          ! Offset to the symmetry block of this shell pair.
          iOff_L = iSSOff(j3,j4)
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Compute index within the irrep. Keep the indexation
          ! of the two basis set sets apart.

          jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)-nBas(j2)
          kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+jOffSO(j3)
          lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+jOffSO(j4)

          nijkl = 0
          do lSOl=lSO,lSO+lBas-1
            iD = iShlSO(lSOl)
            iShD = iSOShl(lSOl)
            nD = nBasSh(j4,iShD)
            do kSOk=kSO,kSO+kBas-1
              iC = iShlSO(kSOk)
              iShC = iSOShl(kSOk)
              nC = nBasSh(j3,iShC)

              if (iShC == iShD) then
                if (j12 == 0) then
                  kl = iTri(iC,iD)
                else if (j3 > j4) then
                  kl = (iC-1)*nD+iD
                else
                  kl = (iD-1)*nC+iC
                end if
              else
                if (iShC >= iShD) then
                  kl = (iD-1)*nC+iC
                else
                  kl = (iC-1)*nD+iD
                end if
              end if

              do jSOj=jSO,jSO+jBas-1
                iAux = jSOj
                nijkl = nijkl+1
                !                                                      *
                !*******************************************************
                !                                                      *
                if ((j12 /= 0) .or. (kSOk >= lSOl)) then
                  kl_B = (iAux-1)*n3C+kl+iOff_L
                  TInt(kl_B) = SOint(nijkl,memSO2)
                end if
                !                                                      *
                !*******************************************************
                !                                                      *
              end do
            end do
          end do

        end do
      end do

    end do
  end do
end do

return

end subroutine IndSft_RI_3
