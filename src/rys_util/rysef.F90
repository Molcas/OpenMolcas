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

implicit real*8(A-H,O-Z)
#include "TriInd.fh"
#include "real.fh"
#include "print.fh"
real*8 xyz2D(nRys,mArg,3,0:neMax,0:nfMax), PreFct(mArg), Scrtch(nRys,mArg), EFInt(nArg,meMin:meMax,mfMin:mfMax)
logical AeqB, CeqD
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
character*80 Label
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
ne = (neMax+1)*(neMax+2)/2
nf = (nfMax+1)*(nfMax+2)/2

if ((ne > IJ_Max) .or. (nf > IJ_Max)) then
  write(6,*) 'ne,nf=',ne,nf
  call WarningMessage(2,'Increase IJ_Max to the larger of the above!')
  call Abend()
end if

do ief=1,ne*nf
  if = (ief-1)/ne+1
  ie = ief-(if-1)*ne

  ixe = iTriInd(1,ie)
  iye = iTriInd(2,ie)
  ixye = ixe+iye

  ixf = iTriInd(1,if)
  iyf = iTriInd(2,if)
  ixyf = ixf+iyf
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nzeMax = max(0,neMax-ixe-iye)
  nzfMax = max(0,nfMax-ixf-iyf)
  nzeMin = max(0,neMin-ixe-iye)
  if (AeqB) nzeMin = nzeMax
  nzfMin = max(0,nfMin-ixf-iyf)
  if (CeqD) nzfMin = nzfMax

  nItem = (nzeMax-nzeMin+1)*(nzfMax-nzfMin+1)
  if (nItem > 1) then

    ! Precompute for all arguments Ix*Iy, avoid multiplying with ones.

    ! Combine with possible Iz

    if (ixe+ixf+iye+iyf == 0) then

      call RysEF1(xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin, &
                  nzeMax,nzfMin,nzfMax)

    else if (ixe+ixf == 0) then

      call RysEF0(xyz2D(1,1,2,iye,iyf),xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf, &
                  ixye,ixyf,nzeMin,nzeMax,nzfMin,nzfMax)

    else if (iye+iyf == 0) then

      call RysEF0(xyz2D(1,1,1,ixe,ixf),xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf, &
                  ixye,ixyf,nzeMin,nzeMax,nzfMin,nzfMax)

    else

      do iArg=1,mArg
        do iRys=1,nRys
          Scrtch(iRys,iArg) = xyz2D(iRys,iArg,1,ixe,ixf)*xyz2D(iRys,iArg,2,iye,iyf)
        end do
      end do
      call RysEF0(Scrtch,xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf, &
                  nzeMin,nzeMax,nzfMin,nzfMax)

    end if

  else

    ! Here if only one triplet of 2D-integrals

    ! Contract over roots

    if (ixe+ixf+iye+iyf == 0) then

      call RysEF2(xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin, &
                  nzeMax,nzfMin,nzfMax)

    else if (ixe+ixf == 0) then

      call RysEF3(xyz2D(1,1,2,iye,iyf),xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf, &
                  ixye,ixyf,nzeMin,nzeMax,nzfMin,nzfMax)

    else if (iye+iyf == 0) then

      call RysEF3(xyz2D(1,1,1,ixe,ixf),xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf, &
                  ixye,ixyf,nzeMin,nzeMax,nzfMin,nzfMax)

    else

      call RysEF4(xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin, &
                  nzeMax,nzfMin,nzfMax)

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
