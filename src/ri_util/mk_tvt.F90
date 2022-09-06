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

subroutine Mk_tVt(TInt,nTheta_All,tVt,nTheta,List2,mData,iPrm,nPrm,iAng,jAng,nk,nl,Indkl,nkl,iAL,nA,nB)

implicit real*8(a-h,o-z)
real*8 TInt(nTheta_All,nTheta_All), tVt(nTheta,nTheta)
integer List2(mData,nTheta_All), iPrm(nPrm), Indkl(nkl), iAL(nA,nB)

!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Mk_tVt: TInt',' ',TInt,nTheta_All,nTheta_All)
call iVcPrt('iPrm',' ',iPrm,nPrm)
call iVcPrt('Indkl',' ',Indkl,nkl)
call iVcPrt('iAL',' ',iAL,nA*nB)
#endif
call FZero(tVt,nTheta**2)
iA = iAng+1
jA = jAng+1
do iTheta_All=1,nTheta_All
  kComp = List2(3,iTheta_All)
  lComp = List2(4,iTheta_All)
  ik = List2(5,iTheta_All)
  il = List2(6,iTheta_All)
  if (iAng == jAng) then
    iTheta_Full = ik*(ik-1)/2+il
  else
    iTheta_Full = (il-1)*nk+ik
  end if
  !if ((iAL(kComp,lComp) == 1) .and. (kComp == iA) .and. (lComp == jA) .and. (iPrm(iTheta_Full) == 1)) then
  if ((kComp == iA) .and. (lComp == jA) .and. (iPrm(iTheta_Full) == 1)) then
    iTheta = Indkl(iTheta_Full)

    do jTheta_All=1,nTheta_All
      mComp = List2(3,jTheta_All)
      nComp = List2(4,jTheta_All)
      im = List2(5,jTheta_All)
      in = List2(6,jTheta_All)
      if (iAng == jAng) then
        jTheta_Full = im*(im-1)/2+in
      else
        jTheta_Full = (in-1)*nk+im
      end if
      !if ((iAL(mComp,nComp) == 1) .and. (mComp == iA) .and. (nComp == jA) .and. (iPrm(jTheta_Full) == 1)) then
      if ((mComp == iA) .and. (nComp == jA) .and. (iPrm(jTheta_Full) == 1)) then
        jTheta = Indkl(jTheta_Full)

        tVt(iTheta,jTheta) = tVt(iTheta,jTheta)+TInt(iTheta_All,jTheta_All)
        !tVt(iTheta,jTheta) = tVt(iTheta,jTheta)+abs(TInt(iTheta_All,jTheta_All))

      end if
    end do

  end if
end do

#ifdef _DEBUGPRINT_
call RecPrt('tVt',' ',tVt,nTheta,nTheta)
#else
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(iAL)
#endif

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nl)

end subroutine Mk_tVt
