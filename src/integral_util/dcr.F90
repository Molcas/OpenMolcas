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

subroutine DCR(Lambda,iStab1,nStab1,iStab2,nStab2,iDCR,mDCR)
!***********************************************************************
! Oject: to compute the double coset representatives (DCR) and Lambda. *
!                                                                      *
! Called from: OneEl                                                   *
!              NucAtt                                                  *
!              TwoEl                                                   *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iOper
use dcr_mod, only: nIndex, Done

implicit none
integer nStab1, nStab2
integer iStab1(0:nStab1-1), iStab2(0:nStab2-1), iDCR(0:7)
integer, save :: index(50), Lambda_all(1275), mDCR_all(1275), iDCR_all(0:7,1275)
integer mDCR
integer Ind1, i, iIrrep, k, Ind2, ij, Lambda

Ind1 = 0
do i=1,nStab1-1
  do iIrrep=1,nIrrep-1
    if (iStab1(i) == iOper(iIrrep)) then
      Ind1 = Ind1+2**(iIrrep-1)
      Go To 11
    end if
  end do
11 continue
end do
do k=1,nIndex
  if (Ind1 == index(k)) then
    Ind1 = k
    Go To 10
  end if
end do
nIndex = nIndex+1
index(nIndex) = Ind1
Ind1 = nIndex
10 continue

Ind2 = 0
do i=1,nStab2-1
  do iIrrep=1,nIrrep-1
    if (iStab2(i) == iOper(iIrrep)) then
      Ind2 = Ind2+2**(iIrrep-1)
      Go To 21
    end if
  end do
21 continue
end do
do k=1,nIndex
  if (Ind2 == index(k)) then
    Ind2 = k
    Go To 20
  end if
end do
nIndex = nIndex+1
index(nIndex) = Ind2
Ind2 = nIndex
20 continue

ij = max(Ind1,Ind2)*(max(Ind1,Ind2)-1)/2+min(Ind1,Ind2)

if (.not. Done(ij)) then
  call DCR_Internal(Lambda_all(ij),iStab1,nStab1,iStab2,nStab2,iDCR_all(0,ij),mDCR_all(ij))
  Done(ij) = .true.
end if
Lambda = Lambda_all(ij)
mDCR = mDCR_all(ij)
call ICopy(mDCR,iDCR_all(0,ij),1,iDCR,1)

contains

subroutine DCR_Internal(Lambda,iStab1,nStab1,iStab2,nStab2,iDCR,mDCR)

  use Symmetry_Info, only: nIrrep, iOper

  implicit none
  integer Lambda, nStab1, nStab2, mDCR
  integer iStab1(0:nStab1-1), iStab2(0:nStab2-1), iDCR(0:7)
  integer iScrt(0:7,0:7)

  integer i, j, k, jik

  iScrt(:,:) = 0

  ! Construct all UGV subgroups. We accumulate the number of times an
  ! operator appears in each UGV subgroup.

  do i=0,nIrrep-1
    do j=0,nStab1-1
      do k=0,nStab2-1
        jik = ieor(iStab1(j),ieor(iOper(i),iStab2(k)))
        iScrt(i,jik) = iScrt(i,jik)+1
      end do
    end do
  end do

  ! Find Lambda. Lambda is the number of times an operator is in the
  ! subgroup UGV. Look only in the first subgroup.

  do i=0,7
    if (iScrt(0,i) /= 0) Lambda = iScrt(0,i)
  end do

  ! Find the unique double cosets (DC) and construct the double coset
  ! representatives (DCR)

  ! Move the first double coset representative, i.e. take the first
  ! operator which appears in the subgroup UGV.

  mDCR = 0
  do i=0,7
    if (iScrt(0,iOper(i)) /= 0) then
      iDCR(mDCR) = iOper(i)
      mDCR = mDCR+1
      Go To 201
    end if
  end do
201 continue

  ! Now look through the remaining subgroups UGV, if any. If a new
  ! unique subgroup is found take a representative from this set.
  ! Observe that the subgroups UGV are either totally identical or
  ! completely disjoint.

  do i=1,nIrrep-1
    do k=0,nIrrep-1
      ! Check if operator exists in UGV.
      if (iScrt(i,iOper(k)) == 0) Go To 211
      do j=0,mDCR-1
        ! See that no element of UGV is already in DCR.
        if (iDCR(j) == iOper(k)) Go To 210
      end do
211   continue
    end do
    ! Here if new unique subgroup UGV was found.
    do k=0,nIrrep-1
      ! Move a representative to the DCR set.
      if (iScrt(i,iOper(k)) /= 0) then
        iDCR(mDCR) = iOper(k)
        mDCR = mDCR+1
        Go To 221
      end if
    end do
221 continue

210 continue
  end do

  return

end subroutine DCR_Internal

end subroutine DCR
