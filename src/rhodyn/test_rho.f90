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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine test_rho (densityt_time, time)
  use rhodyn_data, only: Nstate, threshold, i, j
  use constants, only: auToFs
  implicit none
!
!***********************************************************************
! Purpose: test density matrix on hermicity and positivity (?)
!
!***********************************************************************
!
  complex(8), dimension(:,:) :: densityt_time
  real(8) :: time, abserror

  abserror=0d0
  do i=1,Nstate
    do j=(i+1),Nstate
      if ((abs(dble(densityt_time(i,j))-dble(densityt_time(j,i)))>=&
                      threshold).and.(abs(dble(densityt_time(i,j))-&
                      dble(densityt_time(j,i)))>=abserror)) then
        abserror = abs(dble(densityt_time(i,j))-&
                       dble(densityt_time(j,i)))
      endif
      if ((abs(aimag(densityt_time(i,j))+aimag(densityt_time(j,i)))>=&
          threshold).and.(abs(aimag(densityt_time(i,j))+&
          aimag(densityt_time(j,i)))>=abserror)) then
        abserror = abs(aimag(densityt_time(i,j)) + aimag(densityt_time(j,i)))
      endif
    enddo
  enddo
  if (abserror>=threshold) then
    write(6,'(2(A,X,G28.16,X))')'time=',time*auToFs,'error=',abserror
  endif
end
