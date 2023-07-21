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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_DensDrv(irc,EOcc,EVir,EFro,CMO)
!
! Jonas Bostrom, Feb 2010
!
! Purpose: To Compute MP2 density from Cholesky MO-vectors and
!          decomposed MP2 amplitudes.

use stdalloc

implicit real*8(a-h,o-z)
real*8 EOcc(*), EVir(*), EFro(*), CMO(*)
character*7 ThisNm
character*15 SecNam
parameter(SecNam='ChoMP2g_DensDrv',ThisNm='DensDrv')
integer lWrk
real*8, allocatable :: Wrk(:)

irc = 0

call mma_maxDBLE(lWrk)
! Leave 5% of the memory unallocated
! ----------------------------------
#ifdef _I8_
lWrk = lWrk*19/20
#else
lWrk = lWrk-lWrk/20
#endif
call mma_allocate(Wrk,lWrk,Label='Wrk')
!Wrk(:) = 0.0D0

call ChoMP2g_Reord_R(irc,Wrk,lWrk)

call ChoMP2g_density1(irc,EOcc,EVir,EFro,Wrk,lWrk)
call ChoMP2g_density2(irc,EOcc,EVir,EFro,Wrk,lWrk)

call mma_deallocate(Wrk)

call ChoMP2g_density3(irc,CMO)

end subroutine ChoMP2g_DensDrv
