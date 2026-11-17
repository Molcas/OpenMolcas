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
! Copyright (C) 2018, Quan Phung                                       *
!               2026, Vic Austen                                       *
!***********************************************************************

#include "compiler_features.h"
#ifdef _ENABLE_CHEMPS2_DMRG_

subroutine mkfg3cu4_chemps2(mkF,NLEV,G1,F1,G2,F2,G3,F3,idxG3,nG3,TRANS,chemroot)

use caspt2_module, only: EPSA
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, byte

implicit none
logical(kind=iwp), intent(in) :: mkF, TRANS
integer(kind=iwp), intent(in) :: NLEV, nG3, chemroot
real(kind=wp), intent(in) :: G1(NLEV,NLEV), F1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp), intent(inout) :: F2(NLEV,NLEV,NLEV,NLEV), F3(nG3)
real(kind=wp), intent(out) :: G3(nG3)
integer(kind=byte), intent(in) :: idxG3(6,nG3)
integer(kind=iwp) :: iG3, iw, jt, ju, jv, jx, jy, jz
real(kind=wp) :: EASUM_CHEMPS2
real(kind=wp), allocatable :: G3T(:,:,:,:,:,:)
real(kind=wp), external :: CU4F3H

! The compact 3-RDM is used directly by CASPT2, while the full 3-RDM is
! needed for the G3-dependent terms of the cumulant reconstruction.
call chemps2_load3pdm(NLEV,idxG3,nG3,G3,.true.,EPSA,F2,chemroot,TRANS)
call mma_allocate(G3T,NLEV,NLEV,NLEV,NLEV,NLEV,NLEV,Label='G3T_CheMPS2')
call chemps2_load3pdm_all(NLEV,G3T,chemroot,TRANS)

EASUM_CHEMPS2 = Zero
do iw=1,NLEV
  EASUM_CHEMPS2 = EASUM_CHEMPS2+EPSA(iw)*G1(iw,iw)
end do

if (mkF) then
  do iG3=1,nG3
    jt = idxG3(1,iG3)
    ju = idxG3(2,iG3)
    jv = idxG3(3,iG3)
    jx = idxG3(4,iG3)
    jy = idxG3(5,iG3)
    jz = idxG3(6,iG3)

    F3(iG3) = F3(iG3)+EASUM_CHEMPS2*G3(iG3)
    do iw=1,NLEV
      F3(iG3) = F3(iG3)-Half*EPSA(iw)*(G1(jt,iw)*G3T(iw,jv,jy,ju,jx,jz)+G1(iw,ju)*G3T(jt,jv,jy,iw,jx,jz)+ &
                                       G1(jv,iw)*G3T(iw,jt,jy,jx,ju,jz)+G1(iw,jx)*G3T(jv,jt,jy,iw,ju,jz)+ &
                                       G1(jy,iw)*G3T(iw,jt,jv,jz,ju,jx)+G1(iw,jz)*G3T(jy,jt,jv,iw,ju,jx))
    end do

    F3(iG3) = F3(iG3)+CU4F3H(NLEV,EPSA,EASUM_CHEMPS2,G1,G2,F1,F2,jt,ju,jv,jx,jy,jz)
  end do
end if

call mma_deallocate(G3T)

end subroutine mkfg3cu4_chemps2

#elif ! defined (EMPTY_FILES)
#include "macros.fh"
dummy_empty_procedure(mkfg3cu4_chemps2)
#endif
