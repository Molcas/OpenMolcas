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

subroutine IndSft_Cho_Diag(TInt,lInt,iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOint,nSOint)
!***********************************************************************
!  object: to sift and index the SO integrals.                         *
!                                                                      *
!          the indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          april '90                                                   *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: nIrrep
use Index_Functions, only: iTri, nTri_Elem
use SOAO_Info, only: iAOtSO, iOffSO
use Cholesky, only: iShlSO, iSOShl, nBstSh, ShA, ShB
use sort_data, only: nSkip
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lInt, iCmp(4), iShell(4), iBas, jBas, kBas, lBas, iAO(4), iAOst(4), ijkl, nSOint
real(kind=wp), intent(inout) :: TInt(lInt)
logical(kind=iwp), intent(in) :: Shijij
real(kind=wp), intent(in) :: SOint(ijkl,nSOint)
#include "print.fh"
integer(kind=iwp) :: i1, i12, i2, i3, i34, i4, irout, ISHLI, ISHLJ, iSO, iSOi, iSOij, iSOkl, iSym(0:7), ix, j, j1, j12, j2, j2max, &
                     j3, j4, jCmpMx, jprint, jSO, jSOj, jSym(0:7), k12, k34, KIJ, kSO, kSOk, kSym(0:7), lCmpMx, lSO, lSOl, &
                     lSym(0:7), memSO2, nijkl, NUMI, NUMJ
real(kind=wp) :: r1, r2, tr1 = Zero, tr2 = Zero
logical(kind=iwp) :: qij, qijij, qkl, Shij, Shkl
real(kind=wp), external :: ddot_

irout = 39
jprint = nprint(irout)
k12 = 0
k34 = 0
if (jPrint >= 49) then
  r1 = DDot_(ijkl*nSOInt,SOInt,1,[One],0)
  r2 = DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
  tr1 = tr1+r1
  tr2 = tr2+r2
  write(u6,*) ' Sum=',r1,tr1
  write(u6,*) ' Dot=',r2,tr2
end if
if (jprint >= 99) call RecPrt(' in indsft:SOint ',' ',SOint,ijkl,nSOint)
memSO2 = 0

! allocate space to store integrals to gether with their
! Symmetry batch and sequence number
! To avoid conflicts in using memory this is done in the
! subroutine PSOAO

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.

Shij = iShell(1) == iShell(2)
Shkl = iShell(3) == iShell(4)
do i1=1,iCmp(1)
  do j=0,nIrrep-1
    ix = 0
    if (iAOtSO(iAO(1)+i1,j) > 0) ix = 2**j
    iSym(j) = ix
  end do
  jCmpMx = iCmp(2)
  if (Shij) jCmpMx = i1
  do i2=1,jCmpMx
    do j=0,nIrrep-1
      ix = 0
      if (iAOtSO(iAO(2)+i2,j) > 0) ix = 2**j
      jSym(j) = ix
    end do
    qij = i1 == i2
    if (iShell(2) > iShell(1)) then
      i12 = iCmp(2)*(i1-1)+i2
    else
      i12 = iCmp(1)*(i2-1)+i1
    end if
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
        if (iShell(4) > iShell(3)) then
          i34 = iCmp(4)*(i3-1)+i4
        else
          i34 = iCmp(3)*(i4-1)+i3
        end if
        if (Shijij .and. (i34 > i12)) cycle
        qijij = Shijij .and. (i12 == i34)

        ! loop over Irreps which are spanned by the basis function.
        ! again, the loop structure is restricted to ensure unique
        ! integrals.

        do j1=0,nIrrep-1
          if (iSym(j1) == 0) cycle
          j2max = nIrrep-1
          if (Shij .and. qij) j2max = j1
          do j2=0,j2max
            if (jSym(j2) == 0) cycle
            j12 = ieor(j1,j2)
            if (qijij) then
              if (Shij .and. qij) then
                k12 = nTri_Elem(j1)+j2+1
              else if (Shij) then
                k12 = nIrrep*j1+j2+1
              else if (iShell(1) > iShell(2)) then
                k12 = nIrrep*j1+j2+1
              else
                k12 = nIrrep*j2+j1+1
              end if
            end if

            do j3=0,nIrrep-1
              if (kSym(j3) == 0) cycle
              j4 = ieor(j12,j3)
              if (lSym(j4) == 0) cycle
              if (Shkl .and. qkl .and. (j4 > j3)) cycle
              if (qijij) then
                if (Shkl .and. qkl) then
                  k34 = nTri_Elem(j3)+j4+1
                else if (Shkl) then
                  k34 = nIrrep*j3+j4+1
                else if (iShell(3) > iShell(4)) then
                  k34 = nIrrep*j3+j4+1
                else
                  k34 = nIrrep*j4+j3+1
                end if
                if (k34 > k12) cycle
              end if

              memSO2 = memSO2+1
              if ((nSkip(j1+1)+nSkip(j2+1)+nSkip(j3+1)+nSkip(j4+1)) /= 0) cycle

              ! Compute absolute starting SO index
              iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)+iOffSO(j1)
              jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)+iOffSO(j2)
              kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+iOffSO(j3)
              lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+iOffSO(j4)

              nijkl = 0
              do lSOl=lSO,lSO+lBas-1
                do kSOk=kSO,kSO+kBas-1
                  iSOkl = iTri(kSOk,lSOl)
                  do jSOj=jSO,jSO+jBas-1
                    do iSOi=iSO,iSO+iBas-1
                      nijkl = nijkl+1
                      iSOij = iTri(iSOi,jSOj)

                      if (iSOij == iSOkl) then
                        ISHLI = ISOSHL(ISOI)
                        ISHLJ = ISOSHL(JSOJ)
                        NUMI = NBSTSH(ISHLI)
                        NUMJ = NBSTSH(ISHLJ)
                        if ((ISHLI == ISHLJ) .and. (ISHLI == SHA)) then
                          KIJ = ITRI(ISHLSO(ISOI),ISHLSO(JSOJ))
                        else
                          if ((ISHLI == SHA) .and. (ISHLJ == SHB)) then
                            KIJ = NUMI*(ISHLSO(JSOJ)-1)+ISHLSO(ISOI)
                          else if ((ISHLJ == SHA) .and. (ISHLI == SHB)) then
                            KIJ = NUMJ*(ISHLSO(ISOI)-1)+ISHLSO(JSOJ)
                          else
                            call CHO_QUIT('Integral error',104)
                            KIJ = -999999
                          end if
                        end if
                        TInt(KIJ) = SOint(nijkl,memSO2)
                      end if

                    end do
                  end do
                end do
              end do

            end do
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine IndSft_Cho_Diag
