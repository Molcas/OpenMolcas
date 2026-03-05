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
! x:        global scalar
! op:       global operation '+','*','max','min','absmax','absmin'
subroutine GADGOP_Scal(x,op)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: wp

implicit none
real(kind=wp) :: x
character(len=*) :: op
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"

if (Is_Real_Par()) call ga_dgop(MT_DBL,x,1,op)
#else
#include "macros.fh"
unused_var(x)
unused_var(op)
#endif

end subroutine GADGOP_Scal
