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
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine test_rho(densityt_time,time)
!***********************************************************************
! Purpose: test density matrix on hermicity and positivity (?)
!***********************************************************************

use rhodyn_data, only: Nstate, threshold
use Constants, only: Zero, auToFs
use Definitions, only: wp, iwp, u6

implicit none
complex(kind=wp), intent(in) :: densityt_time(:,:)
real(kind=wp), intent(in) :: time
integer(kind=iwp) :: i, j
real(kind=wp) :: abserror

abserror = Zero
do i=1,Nstate
  do j=(i+1),Nstate
    if ((abs(real(densityt_time(i,j))-real(densityt_time(j,i))) >= threshold) .and. &
        (abs(real(densityt_time(i,j))-real(densityt_time(j,i))) >= abserror)) then
      abserror = abs(real(densityt_time(i,j))-real(densityt_time(j,i)))
    end if
    if ((abs(aimag(densityt_time(i,j))+aimag(densityt_time(j,i))) >= threshold) .and. &
        (abs(aimag(densityt_time(i,j))+aimag(densityt_time(j,i))) >= abserror)) then
      abserror = abs(aimag(densityt_time(i,j))+aimag(densityt_time(j,i)))
    end if
  end do
end do
if (abserror >= threshold) then
  write(u6,'(2(a,1x,g28.16,1x))') 'time=', time*auToFs, 'error=', abserror
end if

end subroutine test_rho
