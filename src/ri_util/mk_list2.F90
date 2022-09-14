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

subroutine Mk_List2(List2,nTheta_All,mData,nSO_Tot,iCnttp,nTest,ijS_req)

use Index_Functions, only: iTri, nTri_Elem1
use Basis_Info, only: dbsc, Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nTheta_All, mData, nSO_Tot, iCnttp, nTest, ijS_req
integer(kind=iwp), intent(out) :: List2(2*mData,nTheta_All)
integer(kind=iwp) :: iAng, iAng_, iCmp, iCmp_, iCont, iCont_, iiSO, ijS, ijSO, iShll, iShll_, iSO, iSO_, jAng, jAng_, jCmp_, &
                     jCont_, jjSO, jShll, jShll_, jSO, jSO_Max, mCmp, mSO, nCmp, nCont, nSO
logical(kind=iwp) :: Only_DB
integer(kind=iwp), allocatable :: iList(:,:)

call mma_allocate(iList,mData,nSO_Tot,Label='iList')

Only_DB = ijS_req /= 0

! Generate intermediate list
! Generate list2, shell blocked!

ijSO = 0
iiSO = 0
iSO_ = 0
do iAng=0,nTest
  iShll = dbsc(iCnttp)%iVal+iAng
  nCmp = nTri_Elem1(iAng)
  if (Shells(iShll)%Prjct) nCmp = 2*iAng+1
  nSO = nCmp*Shells(iShll)%nBasis
  do iCmp=1,nCmp
    nCont = Shells(iShll)%nBasis
    do iCont=1,nCont
      iSO_ = iSO_+1
      iList(1,iSO_) = iAng
      iList(2,iSO_) = iCmp
      iList(3,iSO_) = iCont
      iList(4,iSO_) = iShll
    end do
  end do
  !write(u6,*) 'iSO_=',iSO_

  jjSO = 0
  do jAng=0,iAng
    !write(u6,*) 'iAng,jAng=',iAng,jAng
    jShll = dbsc(iCnttp)%iVal+jAng
    mCmp = nTri_Elem1(jAng)
    if (Shells(jShll)%Prjct) mCmp = 2*jAng+1

    mSO = mCmp*Shells(jShll)%nBasis

    ijS = iTri(iAng+1,jAng+1)

    if ((.not. Only_DB) .or. (ijS == ijS_req)) then
      do iSO=iiSO+1,iiSO+nSO
        iAng_ = iList(1,iSO)
        iCmp_ = iList(2,iSO)
        iCont_ = iList(3,iSO)
        iShll_ = iList(4,iSO)

        jSO_Max = jjSO+mSO
        if (jAng == iAng) jSO_Max = iSO
        do jSO=jjSO+1,jSO_Max
          ijSO = ijSO+1
          jAng_ = iList(1,jSO)
          jCmp_ = iList(2,jSO)
          jCont_ = iList(3,jSO)
          jShll_ = iList(4,jSO)

          !write(u6,*) 'iSO,jSO,ijSO=',iSO,jSO,ijSO
          List2(1,ijSO) = iAng_
          List2(2,ijSO) = jAng_
          List2(3,ijSO) = iCmp_
          List2(4,ijSO) = jCmp_
          List2(5,ijSO) = iCont_
          List2(6,ijSO) = jCont_
          List2(7,ijSO) = iShll_
          List2(8,ijSO) = jShll_
        end do   ! jSO
      end do     ! iSO

    end if

    jjSO = jjSO+mSO
  end do         ! jAng

  iiSO = iiSO+nSO
end do           ! iAng

!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*) 'List2'
write(u6,*) '  iAng,  jAng,  iCmp,  jCmp, iCont, jCont, iShll, jShll'
do ijSO=1,nTheta_All
  write(u6,'(8I7)') (List2(i,ijSO),i=1,8)
end do
#endif

call mma_deallocate(iList)

return

end subroutine Mk_List2
