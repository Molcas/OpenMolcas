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

subroutine complete_ext_loop()
!*****************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: icano_nnend, icano_nnsta, indx, isegdownwei, isegsta, isegupwei, mcroot, value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ilpvalue, irot, irtidx, iupwei, lwei, mm, mm0, mmtmp, nn, nntmp
real(kind=wp) :: vetmp, vlptmp, vlptmp1

!write(u6,*) ' ext_test '
!return
do irot=1,mcroot
  irtidx = indx(irot)
  lwei = isegsta !+irtidx

  do iupwei=1,isegupwei
    ilpvalue = 0
    mm0 = lwei
    nn = lwei+icano_nnsta-1
    do nntmp=icano_nnsta,icano_nnend
      nn = nn+1
      mm = mm0
      vlptmp1 = vector1(nn+irtidx)
      vlptmp = vector2(nn+irtidx)
      do mmtmp=1,nntmp-1
        ilpvalue = ilpvalue+1
        vetmp = value_lpext(ilpvalue)
        mm = mm+1
        vector2(mm+irtidx) = vector2(mm+irtidx)+vlptmp1*vetmp
        vlptmp = vlptmp+vector1(mm+irtidx)*vetmp
      end do
      vector2(nn+irtidx) = vlptmp
    end do
    lwei = lwei+isegdownwei
  end do

end do
!...end of complete_ext_loop

end subroutine complete_ext_loop

subroutine complete_ext_loop_g()

use gugaci_global, only: dm1tmp, icano_nnend, icano_nnsta, index_lpext, index_lpext1, index_lpext2, isegdownwei, isegsta, &
                         isegupwei, value_lpext, value_lpext1, value_lpext2, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ilpvalue, indexlp, indexlp1, indexlp2, iupwei, lwei, mm, mm0, mmtmp, nn, nntmp
real(kind=wp) :: valuelp, valuelp1, valuelp2

lwei = isegsta
do iupwei=1,isegupwei
  ilpvalue = 0
  mm0 = lwei
  nn = lwei+icano_nnsta-1
  do nntmp=icano_nnsta,icano_nnend
    nn = nn+1
    mm = mm0
    do mmtmp=1,nntmp-1
      ilpvalue = ilpvalue+1
      mm = mm+1
      indexlp = index_lpext(ilpvalue)
      if (indexlp /= 0) then
        valuelp = value_lpext(ilpvalue)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      end if
      indexlp1 = index_lpext1(ilpvalue)
      if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ilpvalue)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
      indexlp2 = index_lpext2(ilpvalue)
      if (indexlp2 /= 0) then
        valuelp2 = value_lpext2(ilpvalue)
        dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      end if
    end do
  end do
  lwei = lwei+isegdownwei
end do

end subroutine complete_ext_loop_g
