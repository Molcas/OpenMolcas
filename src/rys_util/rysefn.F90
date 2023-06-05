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
!#define _DEBUGPRINT_

subroutine RysEFn(xyz2D,xyz2Dn,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,AeqB,CeqD,nOrdOp)
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

use Index_Functions, only: iTri_Rev, C3_Ind
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, mArg, nRys, neMin, neMax, nfMin, nfMax, meMin, meMax, mfMin, mfMax, nOrdOp
real(kind=wp), intent(in) :: xyz2D (nRys,mArg,3,  0:neMax+2,0:nfMax+2)
real(kind=wp), intent(in) :: xyz2Dn(nRys,mArg,3,2,0:neMax,       0:nfMax)
real(kind=wp), intent(in) :: PreFct(mArg)
real(kind=wp), intent(out) :: EFInt(nArg,6,meMin:meMax,mfMin:mfMax)
logical(kind=iwp), intent(in) :: AeqB, CeqD
integer(kind=iwp) :: ie, ief, if_, itr(2), ixe, ixf, ixye, ixyf, iye, iyf, ne, nf, nzeMax, nzeMin, nzfMax, nzfMin
integer(kind=iwp) ::  Inde, Indf, ize, izf
integer(kind=iwp) :: iRys
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

  do izf=nzfMin,nzfMax
     Indf = C3_Ind(ixyf+izf,ixf,izf)-1
     do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        EFInt(1:mArg,1,Inde,Indf) = xyz2Dn(1,:,1,2,ixe,ixf)*xyz2D (1,:,2,  iye,iyf)*xyz2D (1,:,3,  ize,izf)
        EFInt(1:mArg,2,Inde,Indf) = xyz2Dn(1,:,1,1,ixe,ixf)*xyz2Dn(1,:,2,1,iye,iyf)*xyz2D (1,:,3,  ize,izf)
        EFInt(1:mArg,3,Inde,Indf) = xyz2Dn(1,:,1,1,ixe,ixf)*xyz2D (1,:,2,  iye,iyf)*xyz2Dn(1,:,3,1,ize,izf)
        EFInt(1:mArg,4,Inde,Indf) = xyz2D (1,:,1,  ixe,ixf)*xyz2Dn(1,:,2,2,iye,iyf)*xyz2D (1,:,3,  ize,izf)
        EFInt(1:mArg,5,Inde,Indf) = xyz2D (1,:,1,  ixe,ixf)*xyz2Dn(1,:,2,1,iye,iyf)*xyz2Dn(1,:,3,1,ize,izf)
        EFInt(1:mArg,6,Inde,Indf) = xyz2D (1,:,1,  ixe,ixf)*xyz2D (1,:,2,  iye,iyf)*xyz2Dn(1,:,3,2,ize,izf)
        do iRys=2,nRys
           EFInt(1:mArg,1,Inde,Indf) = EFInt(1:mArg,1,Inde,Indf) + &
                                       xyz2Dn(1,:,1,2,ixe,ixf)*xyz2D (1,:,2,  iye,iyf)*xyz2D (1,:,3,  ize,izf)
           EFInt(1:mArg,2,Inde,Indf) = EFInt(1:mArg,2,Inde,Indf) + &
                                       xyz2Dn(1,:,1,1,ixe,ixf)*xyz2Dn(1,:,2,1,iye,iyf)*xyz2D (1,:,3,  ize,izf)
           EFInt(1:mArg,3,Inde,Indf) = EFInt(1:mArg,3,Inde,Indf) + &
                                       xyz2Dn(1,:,1,1,ixe,ixf)*xyz2D (1,:,2,  iye,iyf)*xyz2Dn(1,:,3,1,ize,izf)
           EFInt(1:mArg,4,Inde,Indf) = EFInt(1:mArg,4,Inde,Indf) + &
                                       xyz2D (1,:,1,  ixe,ixf)*xyz2Dn(1,:,2,2,iye,iyf)*xyz2D (1,:,3,  ize,izf)
           EFInt(1:mArg,5,Inde,Indf) = EFInt(1:mArg,5,Inde,Indf) + &
                                       xyz2D (1,:,1,  ixe,ixf)*xyz2Dn(1,:,2,1,iye,iyf)*xyz2Dn(1,:,3,1,ize,izf)
           EFInt(1:mArg,6,Inde,Indf) = EFInt(1:mArg,6,Inde,Indf) + &
                                       xyz2D (1,:,1,  ixe,ixf)*xyz2D (1,:,2,  iye,iyf)*xyz2Dn(1,:,3,2,ize,izf)
        end do
        EFInt(1:mArg,1,Inde,Indf) = EFInt(1:mArg,1,Inde,Indf)*PreFct(:)
        EFInt(1:mArg,2,Inde,Indf) = EFInt(1:mArg,2,Inde,Indf)*PreFct(:)
        EFInt(1:mArg,3,Inde,Indf) = EFInt(1:mArg,3,Inde,Indf)*PreFct(:)
        EFInt(1:mArg,4,Inde,Indf) = EFInt(1:mArg,4,Inde,Indf)*PreFct(:)
        EFInt(1:mArg,5,Inde,Indf) = EFInt(1:mArg,5,Inde,Indf)*PreFct(:)
        EFInt(1:mArg,6,Inde,Indf) = EFInt(1:mArg,6,Inde,Indf)*PreFct(:)
     end do
  end do
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
     write(Label,'(A,I3,A,I3,A)') ' In RysEF:xx [',iab,',0|',icd,',0]'
     call RecPrt(Label,' ',EFInt(:,1,iab,icd),1,nArg)
     write(Label,'(A,I3,A,I3,A)') ' In RysEF:xy [',iab,',0|',icd,',0]'
     call RecPrt(Label,' ',EFInt(:,2,iab,icd),1,nArg)
     write(Label,'(A,I3,A,I3,A)') ' In RysEF:xz [',iab,',0|',icd,',0]'
     call RecPrt(Label,' ',EFInt(:,3,iab,icd),1,nArg)
     write(Label,'(A,I3,A,I3,A)') ' In RysEF:yy [',iab,',0|',icd,',0]'
     call RecPrt(Label,' ',EFInt(:,4,iab,icd),1,nArg)
     write(Label,'(A,I3,A,I3,A)') ' In RysEF:yz [',iab,',0|',icd,',0]'
     call RecPrt(Label,' ',EFInt(:,5,iab,icd),1,nArg)
     write(Label,'(A,I3,A,I3,A)') ' In RysEF:zz [',iab,',0|',icd,',0]'
     call RecPrt(Label,' ',EFInt(:,6,iab,icd),1,nArg)
  end do
end do
#endif

return

end subroutine RysEFn
