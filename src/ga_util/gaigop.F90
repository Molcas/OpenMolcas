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

! integer global operation; stub routine to ga_igop...
! k(n):     global vector
! op:       global operation '+','*','max','min','absmax','absmin'
subroutine GAIGOP(k,n,op)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: MT_INT
use Definitions, only: ItoB
use, intrinsic :: iso_c_binding, only: c_int
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(inout) :: k(n)
character(len=*), intent(in) :: op
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iblk
! maximum number of integer values handled by a single GAIGOP (ARMCI) call
! integer(kind=iwp), parameter :: MAXBUF = (huge(1_c_int)-mod(int(huge(1_c_int),kind=iwp),ItoB))/ItoB
integer(kind=iwp), parameter :: MAXBUF = 2**30/ItoB
#endif

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  if (n < MAXBUF) then
    call ga_igop(MT_INT,k,n,op)
  else
    do iblk=0,(n-1)/MAXBUF
      call ga_igop(MT_INT,k(1+MAXBUF*iblk),min(n-MAXBUF*iblk,MAXBUF),op)
    end do
  end if
end if
#else
#include "macros.fh"
unused_var(k)
unused_var(_str(op))
#endif

end subroutine GAIGOP
