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

subroutine Mk_TInt_P(TInt_p,nTInt_p,TP,nTP,List2_p,nList2_p,mData,iAng,jAng,npk,List_TP)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTInt_p, nTP, nList2_p, mData, List2_p(mData,nList2_p), iAng, jAng, npk
real(kind=wp), intent(in) :: TInt_p(nTInt_p,nTInt_p)
real(kind=wp), intent(out) :: TP(nTP,nTP)
integer(kind=iwp), intent(out) :: List_TP(2,nTP)
integer(kind=iwp) :: iA, iList2_p, iTP, jA, jList2_p, jTP, k, kAng, kComp, l, lAng, lComp, m, mAng, mComp, n, nAng, nComp

iA = iAng+1
jA = jAng+1
TP(:,:) = Zero
do iList2_p=1,nList2_p
  kAng = List2_p(1,iList2_p)
  lAng = List2_p(2,iList2_p)
  kComp = List2_p(3,iList2_p)
  lComp = List2_p(4,iList2_p)
  !write(u6,*) 'kComp,lComp=',kComp,lComp
  !if ((kAng == iAng) .and. (iAL(kComp,lComp) == 1) .and. (lAng == jAng) .and. (kComp == iA) .and. (lComp == jA)) then
  if ((kAng == iAng) .and. (lAng == jAng) .and. (kComp == iA) .and. (lComp == jA)) then

    k = List2_p(5,iList2_p)
    l = List2_p(6,iList2_p)
    if (iAng == jAng) then
      iTP = nTri_Elem(k-1)+l
    else
      iTP = (l-1)*npk+k
    end if
    List_TP(1,iTP) = k
    List_TP(2,iTP) = l

    do jList2_p=1,nList2_p
      mAng = List2_p(1,jList2_p)
      nAng = List2_p(2,jList2_p)
      mComp = List2_p(3,jList2_p)
      nComp = List2_p(4,jList2_p)
      !write(u6,*) 'mComp,nComp=',mComp,nComp
      !if ((mAng == iAng) .and. (iAL(mComp,nComp) == 1) .and. (nAng == jAng) .and. (mComp == iA) .and. (nComp == jA)) then
      if ((mAng == iAng) .and. (nAng == jAng) .and. (mComp == iA) .and. (nComp == jA)) then

        m = List2_p(5,jList2_p)
        n = List2_p(6,jList2_p)
        if (iAng == jAng) then
          jTP = nTri_Elem(m-1)+n
        else
          jTP = (n-1)*npk+m
        end if

        TP(iTP,jTP) = TP(iTP,jTP)+TInt_P(iList2_p,jList2_p)
        !TP(iTP,jTP) = TP(iTP,jTP)+abs(TInt_P(iList2_p,jList2_p))

      end if

    end do

  end if
end do

return

end subroutine Mk_TInt_P
