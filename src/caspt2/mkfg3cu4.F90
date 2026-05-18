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
! Copyright (C) 2014, Naoki Nakatani                                   *
!***********************************************************************

#include "compiler_features.h"

#ifdef _ENABLE_BLOCK_DMRG_
subroutine MKFG3CU4(mkF,NLEV,G1,F1,G2,F2,G3,F3,idxG3,nG3,W3)
! Load 1-el, 2-el, and 3-el density matrices to resp. G1, G2, and G3
! and compute 1-el to 4-el contractions of Fock operator F1, F2, and F3.
! For F3, use cumulant reconstruction from 1-el, 2-el, and 3-el density matrices.
!
! Written by N. Nakatani, Oct. 2014

use Symmetry_Info, only: Mul
use constants, only: Half
use definitions, only: iwp, wp, Byte

implicit none
logical(kind=iwp), intent(in) :: mkF,
integer(kind=iwp), intent(in) :: NLEV, nG3
real(kind=wp), intent(out) :: G1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp), intent(out) :: F1(NLEV,NLEV), F2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp), intent(out) :: G3(nG3), F3(nG3)
integer(kind=Byte), intent(in) :: idxG3(6,nG3)
real(kind=wp), intent(in) :: W3(NLEV,NLEV,NLEV,NLEV)
real(kind=wp) G1SUM
integer(kind=iwp) IT, IU, IV, IX, IY, IZ, IW
integer(kind=iwp) JT, JU, JV, JX, JY, JZ
integer(kind=iwp) IZSYM, IYZSYM, IXYZSYM, IVXYZSYM
integer(kind=iwp) IG3
real(kind=wp), external :: CU4F3H

if (NACTEL > 1) then
  ! load 2-el density matrix
  call block_load2pdm(nlev,G2,jstate,jstate)
  ! compute 1-el density matrix from 2-el density matrix
  do iu=1,nlev
    do it=1,nlev
      G1sum = Zero
      if (ism(it) == ism(iu)) then
        do iw=1,nlev
          G1sum = G1sum+G2(iw,iw,it,iu)
        end do
        G1(it,iu) = G1sum/(NACTEL-1)
      end if
    end do
  end do
else
  ! special case for NACTEL = 1
  call block_load1pdm(nlev,G1,jstate,jstate)
end if

do iz=1,nlev
  izSym = ism(iz)
  do iy=1,nlev
    iyzSym = Mul(ism(iy),izSym)
    if (mkF .and. (iyzSym == 1)) then
      do iw=1,nlev
        F1(iy,iz) = F1(iy,iz)+G2(iw,iw,iy,iz)*EPSA(iw)
      end do
    end if
  end do
end do

! skip 3RDM part if NACTEL <= 2
if (NACTEL <= 2) return

do iz=1,nlev
  izSym = ism(iz)
  do iy=1,nlev
    iyzSym = Mul(ism(iy),izSym)
    ! load 3PDM of which is G3(:,:,:,:,iy,iz)
    call block_load3pdm2f(nlev,W3,jstate,jstate,iy,iz)
    if (mkF) then
      do ix=1,nlev
        ixyzSym = Mul(ism(ix),iyzSym)
        do iv=1,nlev
          ivxyzSym = Mul(ism(iv),ixyzSym)
          if (ivxyzSym == 1) then
            do iw=1,nlev
              F2(iv,ix,iy,iz) = F2(iv,ix,iy,iz)+W3(iw,iw,iv,ix)*EPSA(iw)
            end do
          end if
        end do
      end do
    end if

    do iG3=1,NG3
      jt = idxG3(1,iG3)
      ju = idxG3(2,iG3)
      jv = idxG3(3,iG3)
      jx = idxG3(4,iG3)
      jy = idxG3(5,iG3)
      jz = idxG3(6,iG3)
      if ((iy == jy) .and. (iz == jz)) then
        G3(iG3) = W3(jt,ju,jv,jx)
        if (mkF) then
          ! CU4F3 Contrib. :: + G1(lT,lT)*G3(iP,iQ,jP,jQ,kP,kQ)
          F3(iG3) = F3(iG3)+EASUM*G3(iG3)
          do iw=1,nlev
            ! CU4F3 Contrib. :: - 0.5D0*G1(iP,lT)*G3(lT,iQ,jP,jQ,kP,kQ)
            ! CU4F3 Contrib. :: - 0.5D0*G1(lT,iQ)*G3(iP,lT,jP,jQ,kP,kQ)
            F3(iG3) = F3(iG3)-Half*G1(jt,iw)*W3(iw,ju,jv,jx)*EPSA(iw)-Half*G1(iw,ju)*W3(jt,iw,jv,jx)*EPSA(iw)
          end do
        end if
      end if

      if (mkF .and. (iy == jy) .and. (iz == jz)) then
        do iw=1,nlev
          ! CU4F3 Contrib. :: - 0.5D0*G1(jP,lT)*G3(lT,jQ,iP,iQ,kP,kQ)
          ! CU4F3 Contrib. :: - 0.5D0*G1(lT,jQ)*G3(jP,lT,iP,iQ,kP,kQ)
          F3(iG3) = F3(iG3)-Half*G1(jv,iw)*W3(iw,jx,jt,ju)*EPSA(iw)-Half*G1(iw,jx)*W3(jv,iw,jt,ju)*EPSA(iw)
        end do
      end if

      if (mkF .and. (iy == jv) .and. (iz == jx)) then
        do iw=1,nlev
          ! CU4F3 Contrib. :: - 0.5D0*G1(kP,lT)*G3(lT,kQ,iP,iQ,jP,jQ)
          ! CU4F3 Contrib. :: - 0.5D0*G1(lT,kQ)*G3(kP,lT,iP,iQ,jP,jQ)
          F3(iG3) = F3(iG3)-Half*G1(jy,iw)*W3(iw,jz,jt,ju)*EPSA(iw)-Half*G1(iw,jz)*W3(jy,iw,jt,ju)*EPSA(iw)
        end do
      end if
    end do
  end do
end do

if (mkF) then
  do iG3=1,NG3
    it = idxG3(1,iG3)
    iu = idxG3(2,iG3)
    iv = idxG3(3,iG3)
    ix = idxG3(4,iG3)
    iy = idxG3(5,iG3)
    iz = idxG3(6,iG3)
    F3(iG3) = F3(iG3)+CU4F3H(nlev,EPSA,EASUM,G1,G2,F1,F2,it,iu,iv,ix,iy,iz)
  end do
end if

end subroutine MKFG3CU4

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
subroutine empty_MKFG3CU4()
end subroutine empty_MKFG3CU4

#endif
