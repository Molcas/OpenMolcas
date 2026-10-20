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
use Constants, only: Half, Zero
use Definitions, only: wp, iwp, byte
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real(kind=wp) :: EASUM_CHEMPS2
logical(kind=iwp), intent(in) :: mkF, TRANS
integer(kind=iwp), intent(in) :: NLEV, nG3, chemroot
real(kind=wp), intent(in) :: G1(NLEV,NLEV), F1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp), intent(inout) :: F2(NLEV,NLEV,NLEV,NLEV), F3(nG3)
real(kind=wp), intent(out) :: G3(nG3)
integer(kind=byte), intent(in) :: idxG3(6,nG3)
integer(kind=iwp) :: iG3, iw, jt, ju, jv, jx, jy, jz
real(kind=wp), allocatable :: G3T(:)
real(kind=wp), external :: CU4F3H

! The compact 3-RDM is used directly by CASPT2, while the full 3-RDM is
! needed for the G3-dependent terms of the cumulant reconstruction.
call chemps2_load3pdm(NLEV,idxG3,nG3,G3,.true.,EPSA,F2,chemroot,TRANS)
call mma_allocate(G3T,NLEV**6,Label='G3T_CheMPS2')
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
      F3(iG3) = F3(iG3) &
        -Half*G1(jt,iw)*g3val(iw,ju,jv,jx,jy,jz)*EPSA(iw) &
        -Half*G1(iw,ju)*g3val(jt,iw,jv,jx,jy,jz)*EPSA(iw) &
        -Half*G1(jv,iw)*g3val(iw,jx,jt,ju,jy,jz)*EPSA(iw) &
        -Half*G1(iw,jx)*g3val(jv,iw,jt,ju,jy,jz)*EPSA(iw) &
        -Half*G1(jy,iw)*g3val(iw,jz,jt,ju,jv,jx)*EPSA(iw) &
        -Half*G1(iw,jz)*g3val(jy,iw,jt,ju,jv,jx)*EPSA(iw)
    end do

    F3(iG3) = F3(iG3)+CU4F3H(NLEV,EPSA,EASUM_CHEMPS2,G1,G2,F1,F2,jt,ju,jv,jx,jy,jz)
  end do
end if

call mma_deallocate(G3T)

contains

function g3val(p1,q1,p2,q2,p3,q3) result(value)
  integer(kind=iwp), intent(in) :: p1, q1, p2, q2, p3, q3
  integer(kind=iwp) :: idx
  real(kind=wp) :: value

  idx = (p1-1)+NLEV*((p2-1)+NLEV*((p3-1)+NLEV*((q1-1)+NLEV*((q2-1)+NLEV*(q3-1)))))
  value = G3T(1+idx)
end function g3val

end subroutine mkfg3cu4_chemps2

#elif ! defined (EMPTY_FILES)
#include "macros.fh"
dummy_empty_procedure(mkfg3cu4_chemps2)
#endif
