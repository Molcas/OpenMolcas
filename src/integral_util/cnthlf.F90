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
! Copyright (C) 1994, Roland Lindh                                     *
!***********************************************************************

subroutine Cnthlf(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,nPrm2,lZeta,nVec,First,IncVec,A1,A2,A3,Indij)
!***********************************************************************
!                                                                      *
! Object: to do a half transformation. The loop over the two matrix-   *
!         matrix multiplications is segmented such that the end of the *
!         intermediate matrix will not push the start of the same out  *
!         from the cache.                                              *
!                                                                      *
! Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCntr1, nPrm1, nCntr2, nPrm2, lZeta, nVec, IncVec, Indij(lZeta)
real(kind=wp), intent(in) :: Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2), A1(lZeta,nVec)
logical(kind=iwp), intent(in) :: First
real(kind=wp), intent(inout) :: A2(IncVec,nprm2), A3(nVec,nCntr1,nCntr2)
integer(kind=iwp), parameter :: mxnprm = 1000 ! be aware of aCD(fat) basis sets.
integer(kind=iwp) :: ic1, iCntr1, iCntr2, idone(mxnprm), ifirst(mxnprm), iiVec, iPrm1, iPrm2, iZeta, last(mxnprm), minva, mPrm2, &
                     mVec, nnz2(mxnprm), nz2
real(kind=wp) :: C1, C2

if ((nPrm1 > mxnprm) .or. (nPrm2 > mxnprm)) then
  call WarningMessage(2,'CntHlf: nPrm > mxnprm')
  call Abend()
end if

! Sparsity check

nz2 = 0
minva = max(4,(nPrm2+1)/2)
do iCntr2=1,nCntr2
  nnz2(icntr2) = 0
  ifirst(icntr2) = nPrm2+1
  last(icntr2) = 0
  do iPrm2=1,nPrm2
    if (Coeff2(iPrm2,iCntr2) /= Zero) then
      ifirst(icntr2) = min(ifirst(icntr2),iprm2)
      last(icntr2) = max(last(icntr2),iprm2)
      nnz2(icntr2) = nnz2(icntr2)+1
    end if
  end do
  if ((nnz2(icntr2) >= minva) .and. (nz2 == iCntr2-1)) nz2 = iCntr2
end do

! Loop sectioning

do iiVec=1,nVec,IncVec
  mVec = min(IncVec,nVec-iiVec+1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! First quarter transformation

  do iCntr1=1,nCntr1
    idone(1:nprm2) = 0

    do iZeta=1,lZeta
      iPrm2 = (Indij(iZeta)-1)/nPrm1+1
      iPrm1 = Indij(iZeta)-(iPrm2-1)*nPrm1
      C1 = Coeff1(iPrm1,iCntr1)
      if (abs(C1) > Zero) then
        if (idone(iprm2) > 0) then
          A2(1:mVec,iPrm2) = A2(1:mVec,iPrm2)+C1*A1(iZeta,iiVec:iiVec+mVec-1)
        else
          A2(1:mVec,iPrm2) = C1*A1(iZeta,iiVec:iiVec+mVec-1)
          idone(iprm2) = 1
        end if
      end if
    end do ! iZeta
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Second quarter transformation

    do iprm2=1,nprm2
      if (idone(iprm2) == 0) A2(1:mVec,iPrm2) = Zero
    end do

    ic1 = 1
    if (nz2 > 1) then
      if (first) then
        call DGEMM_('n','n',mVec,nz2,nprm2,One,A2,IncVec,Coeff2,nprm2,Zero,A3(iivec,iCntr1,1),nvec*ncntr1)
      else
        call DGEMM_('N','N',mVec,nz2,nPrm2,One,A2,IncVec,Coeff2,nprm2,One,A3(iivec,iCntr1,1),nvec*ncntr1)
      end if
      ic1 = nz2+1
    end if

    do iCntr2=ic1,nCntr2
      if (first) then
        if (nnz2(icntr2) >= minva) then
          call dGeMV_('N',mVec,nPrm2,One,A2,IncVec,Coeff2(1,icntr2),1,Zero,A3(iivec,iCntr1,icntr2),1)
        else
          iprm2 = ifirst(icntr2)
          c2 = coeff2(iprm2,icntr2)
          A3(iiVec:iiVec+mVec-1,iCntr1,iCntr2) = C2*A2(1:mVec,iPrm2)

          iPrm2 = ifirst(icntr2)+1
          mPrm2 = last(iCntr2)-iPrm2+1
          if (mPrm2 > 0) call DNaXpY(mPrm2,mVec,Coeff2(iPrm2,iCntr2),1,A2(1,iPrm2),1,IncVec,A3(iiVec,iCntr1,iCntr2),1,0)

        end if
      else
        if (nnz2(icntr2) >= minva) then
          call dGeMV_('N',mVec,nPrm2,One,A2,IncVec,Coeff2(1,icntr2),1,One,A3(iivec,iCntr1,icntr2),1)
        else
          iPrm2 = ifirst(icntr2)
          mPrm2 = last(icntr2)-iPrm2+1
          call DNaXpY(mPrm2,mVec,Coeff2(iPrm2,iCntr2),1,A2(1,iPrm2),1,IncVec,A3(iiVec,iCntr1,iCntr2),1,0)
        end if
      end if
    end do ! iCntr2
  end do   ! iCntr1

  ! End of loop sectioning

end do     ! iiVec

return

end subroutine Cnthlf
