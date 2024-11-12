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

!#define _DEBUGPRINT_
subroutine IndSft2(iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOint,nSOint,iSOSym,nSOs)
!***********************************************************************
!  object: to sift and index the SO integrals.                         *
!                                                                      *
!          the indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          april '90                                                   *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use k2_arrays, only: Sew_Scr
use SOAO_Info, only: iAOtSO, iOffSO
use lw_Info, only: lwInt, lwSqn, lwSyb
use Gateway_Info, only: ThrInt
use Symmetry_Info, only: nIrrep
use sort_data, only: DimSyB, iStBin, lSll, mxSyP, nSkip, Square, TriSyB
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Constants, only: Zero
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iShell(4), iBas, jBas, kBas, lBas, iAO(4), iAOst(4), ijkl, nSOInt, nSOs, iSOSym(2,nSOs)
logical(kind=iwp), intent(in) :: Shijij
real(kind=wp), intent(in) :: SOint(ijkl,nSOint)
integer(kind=iwp) :: i1, i12, i2, i3, i34, i4, iBin, ij, ipD, ipP1, ipP2, ipP3, ipP4, iqq1, iqq2, iqq3, iqq4, iSO, iSOi, iSq1, &
                     iSq2, iSq3, iSq4, iSqNum, iSyBlk, iSym(0:7), iSym12, iSym34, iSymi, ix, j, j1, j12, j2, j2max, j3, j4, jBin, &
                     jCmpMx, jSO, jSOj, jSqNum, jSyBlk, jSym(0:7), jSymj, k12, k34, kl, kSO, kSOk, kSym(0:7), kSymk, lCmpMx, lSO, &
                     lSOl, lSym(0:7), lSyml, MemSO2, nij, nijkl, nkl, nUt
real(kind=wp) :: A_Int
logical(kind=iwp) :: dupli, qij, qijij, qkl, Shij, Shkl
#ifdef _DEBUGPRINT_
real(kind=wp) :: r1, r2, Tr1 = Zero, Tr2 = Zero
real(kind=wp), external :: DDot_
#endif

k12 = 0
k34 = 0
#ifdef _DEBUGPRINT_
r1 = sum(SOint(:,:))
r2 = DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
tr1 = tr1+r1
tr2 = tr2+r2
write(u6,*) ' Sum=',r1,tr1
write(u6,*) ' Dot=',r2,tr2
call RecPrt(' in indsft:SOint ',' ',SOint,ijkl,nSOint)
#endif
memSO2 = 0

! allocate space to store integrals together with their
! Symmetry batch and sequence number
! To avoid conflicts in using memory this is done in the
! subroutine PSOAO

nUt = -1

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.

Shij = (iShell(1) == iShell(2))
Shkl = (iShell(3) == iShell(4))
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
        qijij = (Shijij .and. (i12 == i34))
        !write(u6,*) 'i1,i2,i3,i4=',i1,i2,i3,i4

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

            iSymi = max(j1,j2)+1
            jSymj = min(j1,j2)+1

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
              !write(u6,*) 'j1,j2,j3,j4=',j1,j2,j3,j4

              memSO2 = memSO2+1
              if ((nSkip(j1+1)+nSkip(j2+1)+nSkip(j3+1)+nSkip(j4+1)) /= 0) cycle

              ! Compute absolute starting SO index
              iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)+iOffSO(j1)
              jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)+iOffSO(j2)
              kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+iOffSO(j3)
              lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+iOffSO(j4)
              !write(u6,*) 'iSO,jSO,kSO,lSO=',iSO,jSO,kSO,lSO

              kSymk = max(j3,j4)+1
              lSyml = min(j3,j4)+1

              iSym12 = TriSyB(iSymi,jSymj)
              iSym34 = TriSyB(kSymk,lSyml)

              ! Order Irrep index canonical
              if (iSym34 > iSym12) then
                iSyBlk = mxSyP*(iSym34-1)+iSym12
                jSyBlk = mxSyP*(iSym12-1)+iSym34
                nij = DimSyB(kSymk,lSyml)
                nkl = DimSyB(iSymi,jSymj)
                iSq1 = nkl
                iSq2 = 1
                iSq3 = 1
                iSq4 = nij
                iQQ1 = lSll(iSyBlk)/nkl
                iQQ2 = nkl+1
                iQQ3 = nij+1
                iQQ4 = lSll(jSyBlk)/nij
                iPP1 = iQQ1
                iPP2 = 0
                iPP3 = 0
                iPP4 = iQQ4
              else
                iSyBlk = mxSyP*(iSym12-1)+iSym34
                jSyBlk = mxSyP*(iSym34-1)+iSym12
                nij = DimSyB(iSymi,jSymj)
                nkl = DimSyB(kSymk,lSyml)
                iSq1 = 1
                iSq2 = nkl
                iSq3 = nij
                iSq4 = 1
                iQQ1 = nkl+1
                iQQ2 = lSll(iSyBlk)/nkl
                iQQ3 = lSll(jSyBlk)/nij
                iQQ4 = nij+1
                iPP1 = 0
                iPP2 = iQQ2
                iPP3 = iQQ3
                iPP4 = 0
              end if
              !write(u6,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
              !write(u6,*) 'iSyBlk,jSyBlk=',iSyBlk,jSyBlk
              !write(u6,*) 'iSq1,iSq2,iSq3,iSq4=',iSq1,iSq2,iSq3,iSq4
              !write(u6,*) 'iQQ1,iQQ2,iQQ3,iQQ4=',iQQ1,iQQ2,iQQ3,iQQ4
              !write(u6,*) 'nkl,nij=',nkl,nij

              ! Duplicate integral if permuted integral is
              ! in the same irrep or if all blocks should
              ! be duplicated.
              dupli = (iSyBlk == jSyBlk) .or. Square

              if (.not. dupli) then

                nijkl = 0
                do lSOl=lSO,lSO+lBas-1
                  do kSOk=kSO,kSO+kBas-1
                    kl = iPD(kSOk,lSOl,iSOSym,nSOs)
                    do jSOj=jSO,jSO+jBas-1
                      do iSOi=iSO,iSO+iBas-1
                        nijkl = nijkl+1
                        A_Int = SOint(nijkl,memSO2)
                        if (abs(A_Int) < ThrInt) cycle
                        ij = iPD(iSOi,jSOj,iSOSym,nSOs)
                        !write(u6,*)
                        !write(u6,*) 'iSOi,jSOj,kSOk,lSOl=',iSOi,jSOj,kSOk,lSO
                        !write(u6,*) 'ij,kl=',ij,kl

                        nUt = nUt+1
                        Sew_Scr(lwInt+nUt) = A_Int
                        iBin = (kl-1)/iQQ1+(ij-1)/iQQ2
                        iSqNum = (kl-iBin*iPP1)*iSq1+(ij-iBin*iPP2)*iSq2-nkl
                        Sew_Scr(lwSqN+nUt) = real(iSqNum,kind=wp)
                        Sew_Scr(lwSyB+nUt) = real(iBin+iStBin(iSyBlk),kind=wp)
                        !write(u6,*) 'iSqNum,iBin=',iSqNum,iBin+iStBin(iSyBlk)

                      end do
                    end do
                  end do
                end do

              else

                nijkl = 0
                do lSOl=lSO,lSO+lBas-1
                  do kSOk=kSO,kSO+kBas-1
                    kl = iPD(kSOk,lSOl,iSOSym,nSOs)
                    do jSOj=jSO,jSO+jBas-1
                      do iSOi=iSO,iSO+iBas-1
                        nijkl = nijkl+1
                        A_Int = SOint(nijkl,memSO2)
                        if (abs(A_Int) < ThrInt) cycle
                        ij = iPD(iSOi,jSOj,iSOSym,nSOs)

                        !write(u6,*)
                        !write(u6,*) 'iSOi,jSOj,kSOk,lSOl=',iSOi,jSOj,kSOk,lSO
                        !write(u6,*) 'ij,kl=',ij,kl

                        nUt = nUt+1
                        Sew_Scr(lwInt+nUt) = A_Int
                        iBin = (kl-1)/iQQ1+(ij-1)/iQQ2
                        iSqNum = (kl-iBin*iPP1)*iSq1+(ij-iBin*iPP2)*iSq2-nkl
                        Sew_Scr(lwSqN+nUt) = real(iSqNum,kind=wp)
                        Sew_Scr(lwSyB+nUt) = real(iBin+iStBin(iSyBlk),kind=wp)
                        !write(u6,*) 'iSqNum,iBin=',iSqNum,iBin+iStBin(iSyBlk)

                        nUt = nUt+1
                        Sew_Scr(lwInt+nUt) = A_Int
                        jBin = (kl-1)/iQQ3+(ij-1)/iQQ4
                        jSqNum = (kl-jBin*iPP3)*iSq3+(ij-jBin*iPP4)*iSq4-nij
                        Sew_Scr(lwSqN+nUt) = real(jSqNum,kind=wp)
                        Sew_Scr(lwSyB+nUt) = real(jBin+iStBin(jSyBlk),kind=wp)
                        !write(u6,*) 'jSqNum,jBin=',jSqNum,jBin+iStBin(jSyBlk)

                      end do
                    end do
                  end do
                end do

              end if

            end do
          end do
        end do

      end do
    end do
  end do
end do

! pass the integral to phase 1 of the bin sorting algorithm

call SORT1A(nUt+1,Sew_Scr(lwInt),Sew_Scr(lwSqN),Sew_Scr(lwSyB))
nUt = 0

return

end subroutine IndSft2
