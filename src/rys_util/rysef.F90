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
! Copyright (C) 1990,1991,1994, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine RysEF(xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,Scrtch,PreFct,AeqB,CeqD)
!***********************************************************************
!                                                                      *
!     Object: to compute integrals corresponding to the primitive set  *
!             used for the HRR. The primitive integrals are generated  *
!             from the 2D-integrals according to the Rys quadrature.   *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90.                                             *
!                                                                      *
!             Modified for kernel routine RysEF0 August '90.           *
!             Modified for kernel routines RysS1, RysS2, and RysS3     *
!             September '90.                                           *
!             Modified for improved vectorization August '91.          *
!             Modified for decreased memory access January '94.        *
!***********************************************************************

use Index_Functions, only: iTri_Rev
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, mArg, nRys, neMin, neMax, nfMin, nfMax, meMin, meMax, mfMin, mfMax
real(kind=wp), intent(in) :: xyz2D(nRys,mArg,3,0:neMax,0:nfMax), PreFct(mArg)
real(kind=wp), intent(out) :: EFInt(nArg,meMin:meMax,mfMin:mfMax)
real(kind=wp), intent(inout) :: Scrtch(nRys,mArg)
logical(kind=iwp), intent(in) :: AeqB, CeqD
integer(kind=iwp) :: ie, ief, if_, itr(2), ixe, ixf, ixye, ixyf, iye, iyf, ne, nf, nItem, nzeMax, nzeMin, nzfMax, nzfMin
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
ne = (neMax+1)*(neMax+2)/2
nf = (nfMax+1)*(nfMax+2)/2

do ief=1,ne*nf
  if_ = (ief-1)/ne+1
  ie = ief-(if_-1)*ne

  itr(:) = iTri_Rev(ie)
  ixye = itr(1)-1
  ixe = itr(2)-1
  iye = ixye-ixe

  itr(:) = iTri_Rev(if_)
  ixyf = itr(1)-1
  ixf = itr(2)-1
  iyf = ixyf-ixf
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nzeMax = max(0,neMax-ixye)
  nzfMax = max(0,nfMax-ixyf)
  nzeMin = max(0,neMin-ixye)
  if (AeqB) nzeMin = nzeMax
  nzfMin = max(0,nfMin-ixyf)
  if (CeqD) nzfMin = nzfMax

  nItem = (nzeMax-nzeMin+1)*(nzfMax-nzfMin+1)
  if (nItem > 1) then

    ! Precompute for all arguments Ix*Iy, avoid multiplying with ones.

    ! Combine with possible Iz

    if (ixye+ixyf == 0) then

      call RysEF1(xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin,nzeMax,nzfMin, &
                  nzfMax)

    else if (ixe+ixf == 0) then

      call RysEF0(xyz2D(:,:,2,iye,iyf),xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf, &
                  nzeMin,nzeMax,nzfMin,nzfMax)

    else if (iye+iyf == 0) then

      call RysEF0(xyz2D(:,:,1,ixe,ixf),xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf, &
                  nzeMin,nzeMax,nzfMin,nzfMax)

    else

      Scrtch(:,:) = xyz2D(:,:,1,ixe,ixf)*xyz2D(:,:,2,iye,iyf)
      call RysEF0(Scrtch,xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin,nzeMax, &
                  nzfMin,nzfMax)

    end if

  else

    ! Here if only one triplet of 2D-integrals

    ! Contract over roots

    if (ixye+ixyf == 0) then

      call RysEF2(xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMax,nzfMax)

    else if (ixe+ixf == 0) then

      call RysEF3(xyz2D(:,:,2,iye,iyf),xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf, &
                  nzeMax,nzfMax)

    else if (iye+iyf == 0) then

      call RysEF3(xyz2D(:,:,1,ixe,ixf),xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf, &
                  nzeMax,nzfMax)

    else

      call RysEF4(xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMax,nzfMax)

    end if

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do
!                                                                      *
!***********************************************************************
!                                                                      *

#ifdef _DEBUGPRINT_
do iab=meMin,meMax
  do icd=mfMin,mfMax
    write(Label,'(A,I3,A,I3,A)') ' In RysEF: [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(1,iab,icd),1,nArg)
  end do
end do
#endif

return

end subroutine RysEF
