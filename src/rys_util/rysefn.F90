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

subroutine RysEFn(xyz2D,xyz2Dn,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,AeqB,CeqD)
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

! The _CHECK_R3_TERM_ option will test the 1/R**3 code by computing a single integral as (x**2+y**2+z**2)/r**3 = 1/r
! That is it should generate the regular ERIs.

use Index_Functions, only: iTri_Rev, C3_Ind
#ifdef _CHECK_R3_TERM_
use Constants, only: Zero
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, mArg, nRys, neMin, neMax, nfMin, nfMax, meMin, meMax, mfMin, mfMax
real(kind=wp), intent(in) :: xyz2D(nRys,mArg,3,0:neMax+2,0:nfMax+2), xyz2Dn(nRys,mArg,3,2,0:neMax+1,0:nfMax+1), PreFct(mArg)
#ifdef _CHECK_R3_TERM_
real(kind=wp), intent(out) :: EFInt(nArg,meMin:meMax,mfMin:mfMax)
#else
real(kind=wp), intent(out) :: EFInt(nArg,6,meMin:meMax,mfMin:mfMax)
#endif
logical(kind=iwp), intent(in) :: AeqB, CeqD
integer(kind=iwp) :: ie, ief, if_, Inde, Indf, iRys, itr(2), ixe, ixf, ixye, ixyf, iye, iyf, ize, izf, ne, nf, nzeMax, nzeMin, &
                     nzfMax, nzfMin
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iab, icd
character(len=80) :: Label
#endif

#ifdef _DEBUGPRINT_
do iab=0,neMax+2
  do icd=0,nfMax+2
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2D(x)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2D(:,:,1,iab,icd),nRys,mArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2D(y)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2D(:,:,2,iab,icd),nRys,mArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2D(z)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2D(:,:,3,iab,icd),nRys,mArg)
  end do
end do
do iab=0,neMax
  do icd=0,nfMax
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2Dn(x)(1)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2Dn(:,:,1,1,iab,icd),nRys,mArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2Dn(y)(1)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2Dn(:,:,2,1,iab,icd),nRys,mArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2Dn(z)(1)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2Dn(:,:,3,1,iab,icd),nRys,mArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2Dn(x)(2)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2Dn(:,:,1,2,iab,icd),nRys,mArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2Dn(y)(2)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2Dn(:,:,2,2,iab,icd),nRys,mArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: xyz2Dn(z)(2)(',iab,',',icd,')'
    call RecPrt(Label,' ',xyz2Dn(:,:,3,2,iab,icd),nRys,mArg)
  end do
end do
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
#     ifdef _CHECK_R3_TERM_
      iRys = 1
      EFInt(1:mArg,Inde,Indf) = Zero
      do iRys=1,nRys
        EFInt(1:mArg,Inde,Indf) = EFInt(1:mArg,Inde,Indf)+ &
                                  xyz2Dn(iRys,:,1,2,ixe,ixf)*xyz2D(iRys,:,2,iye,iyf)*xyz2D(iRys,:,3,ize,izf)+ &
                                  xyz2D(iRys,:,1,ixe,ixf)*xyz2Dn(iRys,:,2,2,iye,iyf)*xyz2D(iRys,:,3,ize,izf)+ &
                                  xyz2D(iRys,:,1,ixe,ixf)*xyz2D(iRys,:,2,iye,iyf)*xyz2Dn(iRys,:,3,2,ize,izf)
      end do
      EFInt(1:mArg,Inde,Indf) = EFInt(1:mArg,Inde,Indf)*PreFct(:)
#     else
      EFInt(1:mArg,1,Inde,Indf) = xyz2Dn(1,:,1,2,ixe,ixf)*xyz2D(1,:,2,iye,iyf)*xyz2D(1,:,3,ize,izf)
      EFInt(1:mArg,2,Inde,Indf) = xyz2Dn(1,:,1,1,ixe,ixf)*xyz2Dn(1,:,2,1,iye,iyf)*xyz2D(1,:,3,ize,izf)
      EFInt(1:mArg,3,Inde,Indf) = xyz2Dn(1,:,1,1,ixe,ixf)*xyz2D(1,:,2,iye,iyf)*xyz2Dn(1,:,3,1,ize,izf)
      EFInt(1:mArg,4,Inde,Indf) = xyz2D(1,:,1,ixe,ixf)*xyz2Dn(1,:,2,2,iye,iyf)*xyz2D(1,:,3,ize,izf)
      EFInt(1:mArg,5,Inde,Indf) = xyz2D(1,:,1,ixe,ixf)*xyz2Dn(1,:,2,1,iye,iyf)*xyz2Dn(1,:,3,1,ize,izf)
      EFInt(1:mArg,6,Inde,Indf) = xyz2D(1,:,1,ixe,ixf)*xyz2D(1,:,2,iye,iyf)*xyz2Dn(1,:,3,2,ize,izf)
      do iRys=2,nRys
        EFInt(1:mArg,1,Inde,Indf) = EFInt(1:mArg,1,Inde,Indf)+ &
                                    xyz2Dn(iRys,:,1,2,ixe,ixf)*xyz2D(iRys,:,2,iye,iyf)*xyz2D(iRys,:,3,ize,izf)
        EFInt(1:mArg,2,Inde,Indf) = EFInt(1:mArg,2,Inde,Indf)+ &
                                    xyz2Dn(iRys,:,1,1,ixe,ixf)*xyz2Dn(iRys,:,2,1,iye,iyf)*xyz2D(iRys,:,3,ize,izf)
        EFInt(1:mArg,3,Inde,Indf) = EFInt(1:mArg,3,Inde,Indf)+ &
                                    xyz2Dn(iRys,:,1,1,ixe,ixf)*xyz2D(iRys,:,2,iye,iyf)*xyz2Dn(iRys,:,3,1,ize,izf)
        EFInt(1:mArg,4,Inde,Indf) = EFInt(1:mArg,4,Inde,Indf)+ &
                                    xyz2D(iRys,:,1,ixe,ixf)*xyz2Dn(iRys,:,2,2,iye,iyf)*xyz2D(iRys,:,3,ize,izf)
        EFInt(1:mArg,5,Inde,Indf) = EFInt(1:mArg,5,Inde,Indf)+ &
                                    xyz2D(iRys,:,1,ixe,ixf)*xyz2Dn(iRys,:,2,1,iye,iyf)*xyz2Dn(iRys,:,3,1,ize,izf)
        EFInt(1:mArg,6,Inde,Indf) = EFInt(1:mArg,6,Inde,Indf)+ &
                                    xyz2D(iRys,:,1,ixe,ixf)*xyz2D(iRys,:,2,iye,iyf)*xyz2Dn(iRys,:,3,2,ize,izf)
      end do
      EFInt(1:mArg,1,Inde,Indf) = EFInt(1:mArg,1,Inde,Indf)*PreFct(:)
      EFInt(1:mArg,2,Inde,Indf) = EFInt(1:mArg,2,Inde,Indf)*PreFct(:)
      EFInt(1:mArg,3,Inde,Indf) = EFInt(1:mArg,3,Inde,Indf)*PreFct(:)
      EFInt(1:mArg,4,Inde,Indf) = EFInt(1:mArg,4,Inde,Indf)*PreFct(:)
      EFInt(1:mArg,5,Inde,Indf) = EFInt(1:mArg,5,Inde,Indf)*PreFct(:)
      EFInt(1:mArg,6,Inde,Indf) = EFInt(1:mArg,6,Inde,Indf)*PreFct(:)
#     endif
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
#   ifdef _CHECK_R3_TERM_
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn: [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(:,iab,icd),1,nArg)
#   else
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn:xx [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(:,1,iab,icd),1,nArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn:xy [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(:,2,iab,icd),1,nArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn:xz [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(:,3,iab,icd),1,nArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn:yy [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(:,4,iab,icd),1,nArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn:yz [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(:,5,iab,icd),1,nArg)
    write(Label,'(A,I3,A,I3,A)') ' In RysEFn:zz [',iab,',0|',icd,',0]'
    call RecPrt(Label,' ',EFInt(:,6,iab,icd),1,nArg)
#   endif
  end do
end do
#endif

return

end subroutine RysEFn
