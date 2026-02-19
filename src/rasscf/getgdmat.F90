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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************

subroutine GetGDMat(GDMat)

use Index_Functions, only: iTri, nTri_Elem
use lucia_data, only: DStmp, Dtmp
use Lucia_Interface, only: Lucia_Util
use rasscf_global, only: iAdr15, lRoots, nAc
use general_data, only: JOBIPH, NCONF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GDMat(nTri_Elem(lRoots),NAC,NAC)
integer(kind=iwp) :: CIDisk1, CIDisk2, iOrb, jRoot, kRoot
real(kind=wp), allocatable :: SDtmp(:), TmpD(:), VecL(:), VecR(:)

call mma_allocate(VecL,NConf,Label='VecL')
call mma_allocate(VecR,NConf,Label='VecR')
call mma_allocate(TmpD,NAC**2,Label='TmpD')
call mma_allocate(SDtmp,NAC**2,Label='SDtmp')
SDtmp(:) = DStmp(:)
DStmp(:) = Zero
TmpD(:) = Dtmp(:)
Dtmp(:) = Zero
CIDisk1 = IADR15(4)
do jRoot=1,lRoots
  call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
  CIDisk2 = IADR15(4)
  do kRoot=1,jRoot
    call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
    !write(u6,*) 'VecL and VecR for states',jRoot,kRoot
    !write(u6,*) (VecL(I),I=0,NConf-1)
    !write(u6,*) (VecR(I),I=0,NConf-1)
    call Lucia_Util('Densi',CI_Vector=VecL(:),RVEC=VECR(:))
    !write(u6,*) 'GDMat for states',jRoot,kRoot
    do IOrb=1,NAC
      GDMat(iTri(jRoot,kRoot),:,IOrb) = Dtmp((IOrb-1)*NAC+1:IOrb*NAC)
      !write(u6,'(10(F8.4,2X))') (GDMat(iTri(jRoot,kRoot),JOrb,IOrb),JOrb=1,NAC)
    end do
  end do
end do
DStmp(:) = SDtmp(:)
DTmp(:) = TmpD(:)
call mma_deallocate(SDtmp)
call mma_deallocate(TmpD)
call mma_deallocate(VecL)
call mma_deallocate(VecR)

end subroutine GetGDMat
