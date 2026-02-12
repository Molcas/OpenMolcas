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

use lucia_data, only: DStmp, Dtmp
use Lucia_Interface, only: Lucia_Util
use rasscf_global, only: lRoots, nAc, iAdr15
use general_data, only: JOBIPH, NCONF
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
! Output
real*8, dimension(lRoots*(lRoots+1)/2,NAC,NAC) :: GDMat
integer CIDisk1, CIDisk2
integer NIJ2, jRoot, kRoot, iOrb, jOrb
real*8, allocatable :: SDtmp(:), TmpD(:)
real*8, allocatable :: VecL(:), VecR(:)

call mma_allocate(VecL,NConf,Label='VecL')
call mma_allocate(VecR,NConf,Label='VecR')
call mma_allocate(TmpD,NAC**2,Label='TmpD')
call mma_allocate(SDtmp,NAC**2,Label='SDtmp')
SDtmp(:) = DStmp(:)
DStmp(:) = 0.0d0
TmpD(:) = Dtmp(:)
Dtmp(:) = 0.0d0
CIDisk1 = IADR15(4)
do jRoot=1,lRoots
  call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
  CIDisk2 = IADR15(4)
  do kRoot=1,jRoot
    call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
    !write(6,*) 'VecL and VecR for states',jRoot,kRoot
    !write(6,*) (VecL(I),I=0,NConf-1)
    !write(6,*) (VecR(I),I=0,NConf-1)
    call Lucia_Util('Densi',CI_Vector=VecL(:),RVEC=VECR(:))
    !write(6,*) 'GDMat for states',jRoot,kRoot
    do IOrb=1,NAC
      do JOrb=1,NAC
        NIJ2 = jRoot*(jRoot-1)/2+kRoot
        GDMat(NIJ2,JOrb,IOrb) = Dtmp(JOrb+(IOrb-1)*NAC)
      end do
      !write(6,'(10(F8.4,2X))') (GDMat(NIJ2,IOrb,JOrb),JOrb=1,NAC)
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
