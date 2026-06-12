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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

module NEVPT2_MOD

use Constants, only: Two, Four, Eight
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nAshT

public :: derex2, derex3, dertG2, dertG3, ex2, ex3, nAshT, tG2, tG3

contains

function tG2(it,iu,iy,iz,G1,G2)

  real(kind=wp) :: tG2
  integer(kind=iwp), intent(in) :: it, iu, iy, iz
  real(kind=wp), intent(in) :: G1(NASHT,NASHT), G2(NASHT,NASHT,NASHT,NASHT)

  ! zero delta
  tG2 = G2(iy,it,iz,iu)
  ! one delta
  if (it == iy) tG2 = tG2-Two*G1(iz,iu)
  if (iu == iz) tG2 = tG2-Two*G1(iy,it)
  if (it == iz) tG2 = tG2+G1(iy,iu)
  if (iu == iy) tG2 = tG2+G1(iz,it)

  ! if (.not. diag) return

  !! two delta
  if ((it == iy) .and. (iu == iz)) tG2 = tG2+Four
  if ((it == iz) .and. (iu == iy)) tG2 = tG2-Two

end function tG2

! subroutine Add_MKBNEVB(ITABS_,IXABS_,IUABS_,IYABS_,IVABS_,IZABS_)
!function tG3(it,iu,iv,ix,iy,iz,G1,G2,G3)
function tG3(it,ix,iu,iy,iv,iz,G1,G2,G3)

  real(kind=wp) :: tG3
  integer(kind=iwp), intent(in) :: it, ix, iu, iy, iv, iz
  real(kind=wp), intent(in) :: G1(NASHT,NASHT), G2(NASHT,NASHT,NASHT,NASHT), G3

  ! \tilde{R}_{abc,def}^{(3)}
  !! zero delta
  tG3 = -G3 !! (iz,iy,ix,iv,iu,it)
  !! one delta
  ! if (it == ix) tG3 = tG3+Two*G2(iz,iy,iv,iu)
  ! if (iu == iy) tG3 = tG3+Two*G2(iz,ix,iv,it)
  ! if (iv == iz) tG3 = tG3+Two*G2(iy,ix,iu,it)
  if (it == ix) tG3 = tG3+Two*G2(iz,iv,iy,iu)
  if (iu == iy) tG3 = tG3+Two*G2(iz,iv,ix,it)
  if (iv == iz) tG3 = tG3+Two*G2(iy,iu,ix,it)

  ! if (it == iy) tG3 = tG3-G2(iz,ix,iv,iu)
  ! if (it == iz) tG3 = tG3-G2(ix,iy,iv,iu)
  ! if (iu == ix) tG3 = tG3-G2(iz,iy,iv,it)
  ! if (iu == iz) tG3 = tG3-G2(iy,ix,iv,it)
  ! if (iv == ix) tG3 = tG3-G2(iy,iz,iu,it)
  ! if (iv == iy) tG3 = tG3-G2(iz,ix,iu,it)
  if (it == iy) tG3 = tG3-G2(iz,iv,ix,iu)
  if (it == iz) tG3 = tG3-G2(ix,iv,iy,iu)
  if (iu == ix) tG3 = tG3-G2(iz,iv,iy,it)
  if (iu == iz) tG3 = tG3-G2(iy,iv,ix,it)
  if (iv == ix) tG3 = tG3-G2(iy,iu,iz,it)
  if (iv == iy) tG3 = tG3-G2(iz,iu,ix,it)

  !! two delta
  !  two same spin
  if ((it == ix) .and. (iu == iy)) tG3 = tG3-Four*G1(iz,iv)
  if ((it == ix) .and. (iv == iz)) tG3 = tG3-Four*G1(iy,iu)
  if ((iu == iy) .and. (iv == iz)) tG3 = tG3-Four*G1(ix,it)
  !  one same spin
  if ((it == ix) .and. (iu == iz)) tG3 = tG3+Two*G1(iy,iv)
  if ((it == ix) .and. (iv == iy)) tG3 = tG3+Two*G1(iz,iu)
  if ((iu == iy) .and. (it == iz)) tG3 = tG3+Two*G1(ix,iv)
  if ((iu == iy) .and. (iv == ix)) tG3 = tG3+Two*G1(iz,it)
  if ((iv == iz) .and. (it == iy)) tG3 = tG3+Two*G1(ix,iu)
  if ((iv == iz) .and. (iu == ix)) tG3 = tG3+Two*G1(iy,it)
  !  one crossing
  if ((it == iy) .and. (iu == ix)) tG3 = tG3+Two*G1(iz,iv)
  if ((it == iz) .and. (iv == ix)) tG3 = tG3+Two*G1(iy,iu)
  if ((iu == iz) .and. (iv == iy)) tG3 = tG3+Two*G1(ix,it)

  if ((it == iy) .and. (iu == iz)) tG3 = tG3-G1(ix,iv)
  if ((it == iy) .and. (iv == ix)) tG3 = tG3-G1(iz,iu)
  if ((it == iz) .and. (iu == ix)) tG3 = tG3-G1(iy,iv)
  if ((it == iz) .and. (iv == iy)) tG3 = tG3-G1(ix,iu)
  if ((iu == ix) .and. (iv == iy)) tG3 = tG3-G1(iz,it)
  if ((iu == iz) .and. (iv == ix)) tG3 = tG3-G1(iy,it)

  !! three delta
  if ((it == ix) .and. (iu == iy) .and. (iv == iz)) tG3 = tG3+Eight

  if ((it == ix) .and. (iu == iz) .and. (iv == iy)) tG3 = tG3-Four
  if ((iu == iy) .and. (it == iz) .and. (iv == ix)) tG3 = tG3-Four
  if ((iv == iz) .and. (it == iy) .and. (iu == ix)) tG3 = tG3-Four

  if ((it == iy) .and. (iu == iz) .and. (iv == ix)) tG3 = tG3+Two
  if ((it == iz) .and. (iu == ix) .and. (iv == iy)) tG3 = tG3+Two

end function tG3

function ex2(it,iu,iy,iz,G1,G2)

  real(kind=wp) :: ex2
  integer(kind=iwp), intent(in) :: it, iu, iy, iz
  real(kind=wp), intent(in) :: G1(NASHT,NASHT), G2(NASHT,NASHT,NASHT,NASHT)

  ex2 = G2(it,iu,iy,iz)
  if (iu == iy) ex2 = ex2+G1(it,iz)

end function ex2

function ex3(it,iu,iv,ix,iy,iz,G1,G2,G3)

  use Definitions, only: iwp, wp

  real(kind=wp) :: ex3
  integer(kind=iwp), intent(in) :: it, iu, iv, ix, iy, iz
  real(kind=wp), intent(in) :: G1(NASHT,NASHT), G2(NASHT,NASHT,NASHT,NASHT), G3

  ex3 = G3
  if (iu == iv) ex3 = ex3+G2(it,ix,iy,iz)
  if (ix == iy) ex3 = ex3+G2(it,iu,iv,iz)
  if (iu == iy) ex3 = ex3+G2(it,iz,iv,ix)
  if ((iu == iv) .and. (ix == iy)) ex3 = ex3+G1(it,iz)

end function ex3

subroutine dertG2(it,iu,iy,iz,DG1,DG2,val)

  integer(kind=iwp), intent(in) :: it, iu, iy, iz
  real(kind=wp), intent(inout) :: DG1(NASHT,NASHT), DG2(NASHT,NASHT,NASHT,NASHT)
  real(kind=wp), intent(in) :: val

  ! zero delta
  DG2(iy,it,iz,iu) = DG2(iy,it,iz,iu)+val
  ! one delta
  if (it == iy) DG1(iz,iu) = DG1(iz,iu)-Two*val
  if (iu == iz) DG1(iy,it) = DG1(iy,it)-Two*val
  if (it == iz) DG1(iy,iu) = DG1(iy,iu)+val
  if (iu == iy) DG1(iz,it) = DG1(iz,it)+val

end subroutine dertG2

subroutine dertG3(it,ix,iu,iy,iv,iz,DG1,DG2,DG3,val)

  integer(kind=iwp), intent(in) :: it, ix, iu, iy, iv, iz
  real(kind=wp), intent(inout) :: DG1(NASHT,NASHT), DG2(NASHT,NASHT,NASHT,NASHT), DG3
  real(kind=wp), intent(in) :: val

  ! \tilde{R}_{abc,def}^{(3)}
  !! zero delta
  DG3 = DG3-val !! (iz,iy,ix,iv,iu,it)
  !! one delta
  if (it == ix) DG2(iz,iv,iy,iu) = DG2(iz,iv,iy,iu)+Two*val
  if (iu == iy) DG2(iz,iv,ix,it) = DG2(iz,iv,ix,it)+Two*val
  if (iv == iz) DG2(iy,iu,ix,it) = DG2(iy,iu,ix,it)+Two*val

  if (it == iy) DG2(iz,iv,ix,iu) = DG2(iz,iv,ix,iu)-val
  if (it == iz) DG2(ix,iv,iy,iu) = DG2(ix,iv,iy,iu)-val
  if (iu == ix) DG2(iz,iv,iy,it) = DG2(iz,iv,iy,it)-val
  if (iu == iz) DG2(iy,iv,ix,it) = DG2(iy,iv,ix,it)-val
  if (iv == ix) DG2(iy,iu,iz,it) = DG2(iy,iu,iz,it)-val
  if (iv == iy) DG2(iz,iu,ix,it) = DG2(iz,iu,ix,it)-val

  !! two delta
  !  two same spin
  if ((it == ix) .and. (iu == iy)) DG1(iz,iv) = DG1(iz,iv)-Four*val
  if ((it == ix) .and. (iv == iz)) DG1(iy,iu) = DG1(iy,iu)-Four*val
  if ((iu == iy) .and. (iv == iz)) DG1(ix,it) = DG1(ix,it)-Four*val
  !  one same spin
  if ((it == ix) .and. (iu == iz)) DG1(iy,iv) = DG1(iy,iv)+Two*val
  if ((it == ix) .and. (iv == iy)) DG1(iz,iu) = DG1(iz,iu)+Two*val
  if ((iu == iy) .and. (it == iz)) DG1(ix,iv) = DG1(ix,iv)+Two*val
  if ((iu == iy) .and. (iv == ix)) DG1(iz,it) = DG1(iz,it)+Two*val
  if ((iv == iz) .and. (it == iy)) DG1(ix,iu) = DG1(ix,iu)+Two*val
  if ((iv == iz) .and. (iu == ix)) DG1(iy,it) = DG1(iy,it)+Two*val
  !  one crossing
  if ((it == iy) .and. (iu == ix)) DG1(iz,iv) = DG1(iz,iv)+Two*val
  if ((it == iz) .and. (iv == ix)) DG1(iy,iu) = DG1(iy,iu)+Two*val
  if ((iu == iz) .and. (iv == iy)) DG1(ix,it) = DG1(ix,it)+Two*val

  if ((it == iy) .and. (iu == iz)) DG1(ix,iv) = DG1(ix,iv)-val
  if ((it == iy) .and. (iv == ix)) DG1(iz,iu) = DG1(iz,iu)-val
  if ((it == iz) .and. (iu == ix)) DG1(iy,iv) = DG1(iy,iv)-val
  if ((it == iz) .and. (iv == iy)) DG1(ix,iu) = DG1(ix,iu)-val
  if ((iu == ix) .and. (iv == iy)) DG1(iz,it) = DG1(iz,it)-val
  if ((iu == iz) .and. (iv == ix)) DG1(iy,it) = DG1(iy,it)-val

end subroutine dertG3

subroutine derex2(it,iu,iy,iz,DG1,DG2,val)

  integer(kind=iwp), intent(in) :: it, iu, iy, iz
  real(kind=wp), intent(inout) :: DG1(NASHT,NASHT), DG2(NASHT,NASHT,NASHT,NASHT)
  real(kind=wp), intent(in) :: val

  DG2(it,iu,iy,iz) = DG2(it,iu,iy,iz)+val
  if (iu == iy) DG1(it,iz) = DG1(it,iz)+val

end subroutine derex2

subroutine derex3(it,iu,iv,ix,iy,iz,DG1,DG2,DG3,val)

  integer(kind=iwp), intent(in) :: it, iu, iv, ix, iy, iz
  real(kind=wp), intent(inout) :: DG1(NASHT,NASHT), DG2(NASHT,NASHT,NASHT,NASHT), DG3
  real(kind=wp), intent(in) :: val

  DG3 = DG3+val
  if (iu == iv) DG2(it,ix,iy,iz) = DG2(it,ix,iy,iz)+val
  if (ix == iy) DG2(it,iu,iv,iz) = DG2(it,iu,iv,iz)+val
  if (iu == iy) DG2(it,iz,iv,ix) = DG2(it,iz,iv,ix)+val
  if ((iu == iv) .and. (ix == iy)) DG1(it,iz) = DG1(it,iz)+val

end subroutine derex3

end module NEVPT2_MOD
