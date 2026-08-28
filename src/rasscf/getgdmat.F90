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
use ci_interfaces, only: Mk_T1DM
use rasscf_global, only: iAdr15, lRoots, nAc
use rasscf_files, only: JOBIPH
use general_data, only: NCONF
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: GDMat(nTri_Elem(lRoots),NAC,NAC)
integer(kind=iwp) :: CIDisk1, CIDisk2, iOrb, jRoot, kRoot
real(kind=wp), allocatable :: TmpD(:), VecL(:), VecR(:)

call mma_allocate(VecL,NConf,Label='VecL')
call mma_allocate(VecR,NConf,Label='VecR')
call mma_allocate(TmpD,NAC**2,Label='TmpD')

CIDisk1 = IADR15(4)
do jRoot=1,lRoots
  call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
  CIDisk2 = IADR15(4)
  do kRoot=1,jRoot
    call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
    call Mk_T1DM(VECR,VECL,nConf,TMPD,NAC**2)
    do IOrb=1,NAC
      GDMat(iTri(jRoot,kRoot),:,IOrb) = TmpD((IOrb-1)*NAC+1:IOrb*NAC)
    end do
  end do
end do

call mma_deallocate(TmpD)
call mma_deallocate(VecL)
call mma_deallocate(VecR)

end subroutine GetGDMat
