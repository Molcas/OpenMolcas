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
! Copyright (C) 2016, Sebastian Wouters                                *
!               2016, Quan Phung                                       *
!***********************************************************************

#include "compiler_features.h"

#ifdef _ENABLE_CHEMPS2_DMRG_
subroutine mkfg3chemps2(mkF,NLEV,G1,F1,G2,F2,G3,F3,idxG3,NG3)

use Symmetry_Info, only: Mul
use sguga, only: SGS
use caspt2_module, only: jState, nActel, EPSA, mState
use definitions, only: iwp, wp, Byte, u6

implicit none
logical(kind=iwp), intent(in) :: mkF
integer(kind=iwp), intent(in) :: NLEV, NG3
real(kind=wp), intent(out) :: G1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp), intent(out) :: F1(NLEV,NLEV), F2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp), intent(out) :: G3(NG3), F3(nG3)
integer(kind=Byte), intent(in) :: idxG3(6,nG3)
integer(kind=iwp) IY, IZ, IW
integer(kind=iwp) IYSYM, IXYSYM
integer(kind=iwp) NAC4

if (NACTEL > 1) then
  NAC4 = NLEV*NLEV*NLEV*NLEV
  call chemps2_load2pdm(nlev,G2,MSTATE(JSTATE))
  call two2onerdm(nlev,NACTEL,G2,G1)
else
  write(u6,*) 'FATAL ERROR: DMRG-CASPT2 with CHEMPS2 does not work with NACTEL=1'
end if

! Double checked with CheMPS2::CASPT2::create_f_dots()
do iz=1,nlev
  iySym = SGS%ism(iz)
  do iy=1,nlev
    ixySym = Mul(SGS%ism(iy),iySym)
    if (mkF .and. (ixySym == 1)) then
      F1(iy,iz) = 0.0
      do iw=1,nlev
        F1(iy,iz) = F1(iy,iz)+G2(iw,iw,iy,iz)*EPSA(iw)
      end do
    end if
  end do
end do

if (NACTEL >= 3) then

  if (mkF) call chemps2_load3pdm(nlev,idxG3,NG3,F3,.false.,EPSA,F2,MSTATE(JSTATE))

  call chemps2_load3pdm(nlev,idxG3,NG3,G3,mkF,EPSA,F2,MSTATE(JSTATE))

end if

end subroutine mkfg3chemps2

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
subroutine empty_mkfg3chemps2()
end subroutine empty_mkfg3chemps2

#endif
