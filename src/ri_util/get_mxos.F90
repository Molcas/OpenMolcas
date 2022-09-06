!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Get_mXOs(kOrb,XO,locc,nSkal,nIrrep,nOcc)

use ChoArr, only: nBasSh
use ExTerm, only: CMOi

implicit real*8(a-h,o-z)
integer kOrb, nOcc(nIrrep), nSkal
real*8 XO(locc,nSkal,nIrrep)
#include "cholesky.fh"
#include "choorb.fh"
#include "exterm.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
call FZero(XO,locc*nSkal*nIrrep)

! Loop over irreps

do ir=1,nIrrep

  ! The next block of X_i,mu

  !call RecPrt('X_i,mu',' ',CMOi(kOrb)%SB(ir)%A2,nOcc(ir),nBas(ir))

  ! Loop over all valence shells

  iOff = 0
  do isk=1,nSkal

    ! Loop over all basis functions of this shell in this irrep.

    !write (6,*) 'isk,nBasSh(ir,isk)=',isk,nBasSh(ir,isk)

    do ib=1,nBasSh(ir,isk)
      kb = iOff+ib ! relative SO index in this irrepp

      ! Loop over all the occupied MOs and pick up the largest coefficient for shell isk

      do iok=1,nOcc(ir)
        !write(6,*) 'iok,kb=',iok,kb
        !write(6,*) 'CMOi(kOrb)%SB(ir)%A2(iok,kb)',CMOi(kOrb)%SB(ir)%A2(iok,kb)
        XO(iok,isk,ir) = max(XO(iok,isk,ir),abs(CMOi(kOrb)%SB(ir)%A2(iok,kb)))
      end do
    end do
    iOff = iOff+nBasSh(ir,isk)
  end do
  !call RecPrt('XO(*,*,ir)',' ',XO(1,1,ir),locc,nskal)
end do

return

end subroutine Get_mXOs
