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

subroutine set_nm(exch,ncut,encut_definition,nk,mg,nTempMagn,hmax,w,encut_rate,TempMagn,nM,EM,dbg)

use Constants, only: Zero, cm_s, hPlanck, kBoltzmann, mBohr
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: exch, ncut, encut_definition, nk, mg, nTempMagn
real(kind=wp), intent(in) :: hmax, W(exch), encut_rate, TempMagn(nTempMagn)
integer(kind=iwp), intent(out) :: nM
real(kind=wp), intent(out) :: EM
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i
real(kind=wp) :: diff, T_High
real(kind=wp), parameter :: boltz_k = kBoltzmann/(cm_s*hPlanck), & ! in cm-1*K-1
                            mu_bohr = mBohr/(cm_s*hPlanck) ! in cm-1*T-1

#include "warnings.h"

nM = 1
EM = Zero
diff = Zero
T_High = Zero
if (nTempMagn > 0) T_High = maxval(TempMagn(:))

if (dbg) then
  write(u6,*) 'exch             = ',exch
  write(u6,*) 'ncut             = ',ncut
  write(u6,*) 'encut_definition = ',encut_definition
  write(u6,*) 'nk               = ',nk
  write(u6,*) 'mg               = ',mg
  write(u6,*) 'nM               = ',nM
  write(u6,*) 'nTempMagn        = ',nTempMagn
  write(u6,*) 'hmax             = ',hmax
  write(u6,*) 'encut_rate       = ',encut_rate
  write(u6,*) 'EM               = ',EM
  write(u6,*) 'TempMagn()       = ',TempMagn(:)
  write(u6,*) 'W()              = ',W(:)
end if

if (encut_definition == 1) then

  if (ncut > exch) then
    nm = exch
    em = w(exch)
  else
    nm = ncut
    em = w(nm)
  end if

else if (encut_definition == 2) then

  nm = exch
  em = nk*boltz_k*T_High+mg*mu_bohr*abs(hmax)

  do i=1,exch
    if (i > 1) diff = w(i)-w(i-1)
    if ((w(i) > em) .and. (diff > 1.0e-4_wp)) then
      nm = i-1
      exit
    end if
  end do

else if (encut_definition == 3) then

  nm = exch
  em = w(exch)*encut_rate

  do i=1,exch
    if (i > 1) diff = w(i)-w(i-1)
    if ((w(i) > em) .and. (diff > 1.0e-4_wp)) then
      nm = i-1
      exit
    end if
  end do

else

  write(u6,'(A)') 'something is wrong with "encut_definition"'
  call quit(_RC_INPUT_ERROR_)

end if !encut_definition

return

end subroutine set_nm
