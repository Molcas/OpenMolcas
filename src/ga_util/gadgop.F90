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
! Copyright (C) 1995, Martin Schuetz                                   *
!               1998, Roland Lindh                                     *
!               2000-2015, Steven Vancoillie                           *
!***********************************************************************

! double global operation; stub routine to ga_dgop...
! x(n):     global vector
! op:       global operation '+','*','max','min','absmax','absmin'
subroutine GADGOP(x,n,op)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: MT_DBL
use Definitions, only: RtoB
use, intrinsic :: iso_c_binding, only: c_int
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: x(n)
character(len=*), intent(in) :: op

#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iblk
! maximum number of real values handled by a single GADGOP (ARMCI) call
! compilers complain about huge(xxx)/RtoB with -Werror=integer-division
integer(kind=iwp), parameter :: MAXBUF = (huge(1_c_int)-mod(int(huge(1_c_int),kind=iwp),RtoB))/RtoB
#endif

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  if (n < MAXBUF) then
    call ga_dgop(MT_DBL,x,n,op)
  else
    do iblk=0,(n-1)/MAXBUF
      call ga_dgop(MT_DBL,x(1+MAXBUF*iblk),min(n-MAXBUF*iblk,MAXBUF),op)
    end do
  end if
end if
#else
#include "macros.fh"
unused_var(x)
unused_var(_str(op))
#endif

end subroutine GADGOP
