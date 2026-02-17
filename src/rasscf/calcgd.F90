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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 11, 2022, created this file.               *
!*****************************************************************

! Subroutine relating to generalized 1-e density matrix (GD) called
! in CMSNewton
!     calculating GD with lucia.
subroutine CalcGD(GD,nGD)

use lucia_data, only: Dtmp, DStmp
use Lucia_Interface, only: Lucia_Util
use rasscf_global, only: iADR15, lRoots, NAC
use general_data, only: JOBIPH, NCONF
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nGD
real(kind=wp) :: GD(nGD)
integer(kind=iwp) :: CIDisk1, CIDisk2, IOffNIJ1, IOffNIJ2, ipq, iqp, jRoot, kRoot, NAC2, p, q
real(kind=wp), allocatable :: SDtmp(:), TmpD(:)
real(kind=wp), allocatable, target :: VecL(:), VecR(:)

NAC2 = NAC**2
call mma_allocate(VecL,NConf,Label='VecL')
call mma_allocate(VecR,NConf,Label='VecR')
call mma_allocate(TmpD,NAC**2,Label='TmpD')
call mma_allocate(SDtmp,NAC**2,Label='SDtmp')
SDtmp(:) = DStmp(:)
TmpD(:) = DTmp(:)
CIDisk1 = IADR15(4)
do jRoot=1,lRoots
  call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
  CIDisk2 = IADR15(4)
  do kRoot=1,jRoot-1
    call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
    call Lucia_Util('Densi',CI_Vector=VecL(:),RVec=VecR(:))
    IOffNIJ1 = (lRoots*(jRoot-1)+kRoot-1)*NAC2
    IOffNIJ2 = (lRoots*(kRoot-1)+jRoot-1)*NAC2
    !write(u6,*) 'GD matrix',jRoot,kRoot
    !call RecPrt(' ',' ',Dtmp,NAC,NAC)
    GD(IOffNIJ1+1:IOffNIJ1+NAC2) = Dtmp(1:NAC2)
    do q=1,NAC
      do p=1,NAC
        ipq = (q-1)*NAC+p
        iqp = (p-1)*NAC+q
        GD(IOffNIJ2+iqp) = Dtmp(ipq)
        !GDMat(NIJ2,q,p) = Dtmp(q+(p-1)*NAC)
      end do
    end do
  end do
  kRoot = jRoot
  call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
  call Lucia_Util('Densi',CI_Vector=VecL(:),RVec=VecR(:))
  IOffNIJ1 = (lRoots+1)*(jRoot-1)*NAC2
  !write(u6,*) 'GD matrix',jRoot,kRoot
  !call RecPrt(' ',' ',Dtmp,NAC,NAC)
  GD(IOffNIJ1+1:IOffNIJ1+NAC2) = Dtmp(1:NAC2)
end do
DStmp(:) = SDtmp(:)
Dtmp(:) = TmpD(:)
call mma_deallocate(SDtmp)
call mma_deallocate(TmpD)
call mma_deallocate(VecL)
call mma_deallocate(VecR)

end subroutine CalcGD
