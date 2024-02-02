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

implicit none
#include "warnings.h"
! input data:
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: exch, ncut, encut_definition, nk, mg, nTempMagn
real(kind=8), intent(in) :: hmax, W(exch), encut_rate, TempMagn(nTempMagn)
logical, intent(in) :: dbg
! output data:
integer, intent(out) :: nM
real(kind=8), intent(out) :: EM
! local variables:
integer :: i
real(kind=8) :: diff, T_High
real(kind=8) :: boltz_k, mu_bohr

! Constants:
boltz_k = 0.6950356_wp   ! in cm^-1*K-1
mu_bohr = 0.466864374_wp ! in cm-1*T-1

nM = 1
EM = 0.0_wp
diff = 0.0_wp
T_High = 0.0_wp
if (nTempMagn > 0) T_High = maxval(TempMagn(1:nTempMagn))

if (dbg) then
  write(6,*) 'exch             = ',exch
  write(6,*) 'ncut             = ',ncut
  write(6,*) 'encut_definition = ',encut_definition
  write(6,*) 'nk               = ',nk
  write(6,*) 'mg               = ',mg
  write(6,*) 'nM               = ',nM
  write(6,*) 'nTempMagn        = ',nTempMagn
  write(6,*) 'hmax             = ',hmax
  write(6,*) 'encut_rate       = ',encut_rate
  write(6,*) 'EM               = ',EM
  write(6,*) 'TempMagn()       = ',TempMagn(1:nTempMagn)
  write(6,*) 'W()              = ',W(1:exch)
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
    if ((w(i) > em) .and. (diff > 1.0d-4)) then
      nm = i-1
      go to 309
    end if
  end do

else if (encut_definition == 3) then

  nm = exch
  em = w(exch)*encut_rate

  do i=1,exch
    if (i > 1) diff = w(i)-w(i-1)
    if ((w(i) > em) .and. (diff > 1.0d-4)) then
      nm = i-1
      go to 309
    end if
  end do

else

  write(6,'(A)') 'something is wrong with "encut_definition"'
  call quit(_RC_INPUT_ERROR_)

end if !encut_definition

309 continue

return

end subroutine set_nm
