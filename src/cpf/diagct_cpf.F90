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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine DIAGCT_CPF()

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use cpf_global, only: ICASE, IFIRST, ILIM, IROW, ISAB, ISMAX, JBUF, JSC, JSY, KBUF, KBUFF1, LBUF, MAX11, NCONF, NORBT, NOV, NOV1, &
                      NTIBUF, NVT5
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp) :: JBUF1, KBUF1, LBUF1, NINTGR
integer(kind=iwp), allocatable :: IBUFL(:), INDCAT(:)
real(kind=wp), allocatable :: A1(:), A2(:), FC(:), FIJ(:), FJI(:), TIBUF(:)
real(kind=wp), allocatable, target :: BUF(:), BUFOUT(:)
integer(kind=iwp), pointer :: iBUFBI(:), INDOUT(:)

NCONF = JSC(ILIM)
KBUF1 = ((RtoI+1)*KBUF+2+(RtoI-1))/RtoI
JBUF1 = ((RtoI+1)*JBUF+2+(RtoI-1))/RtoI
LBUF1 = ((RtoI+1)*LBUF+2+(RtoI-1))/RtoI
call mma_allocate(TIBUF,NTIBUF,label='TIBUF')
call mma_allocate(BUFOUT,max(NVT5*JBUF1,NOV*KBUF1,NOV1*LBUF1,MAX11),label='BUFOUT')
call mma_allocate(INDCAT,max(NVT5,NOV,NOV1),label='INDCAT')
call mma_allocate(IBUFL,max(NVT5,NOV,NOV1,25000),label='IBUFL')
call c_f_pointer(c_loc(BUFOUT),INDOUT,[1])
! Initialize sorting buffer, so that automatic detection of
! uninitialized variables does not give false alarms.
BUFOUT(1:NOV*KBUF1) = Zero
! Similar before SORTB_CPF and SORT_CPF.
call mma_allocate(BUF,KBUFF1+2,label='BUF')
call mma_allocate(A1,ISMAX,label='A1')
call mma_allocate(A2,ISMAX,label='A2')
call c_f_pointer(c_loc(BUF),iBUFBI,[1])
call SORTA_CPF(BUFOUT,INDOUT,INDCAT,IBUFL,TIBUF,ISAB,BUF,iBUFBI,A1,A2,NINTGR)
nullify(iBUFBI)
if (IFIRST == 0) then
  BUFOUT(1:NVT5*JBUF) = Zero
  call SORTB_CPF(BUFOUT,INDOUT,INDCAT,IBUFL,TIBUF,A1,A2,ISAB,BUF)
end if
call mma_deallocate(BUF)
call mma_deallocate(A1)
call mma_deallocate(A2)
BUFOUT(1:NOV1*LBUF) = Zero
call mma_allocate(FC,IROW(NORBT+1),label='FC')
call mma_allocate(FIJ,IROW(NORBT+1),label='FIJ')
call mma_allocate(FJI,IROW(NORBT+1),label='FJI')
call SORT_CPF(BUFOUT,INDOUT,INDCAT,IBUFL,FC,FIJ,FJI,TIBUF)
call DIAG_CPF(ICASE,JSY,BUFOUT,FC,FIJ,FJI)
call mma_deallocate(FC)
call mma_deallocate(FIJ)
call mma_deallocate(FJI)
call mma_deallocate(TIBUF)
nullify(INDOUT)
call mma_deallocate(BUFOUT)
call mma_deallocate(INDCAT)
call mma_deallocate(IBUFL)

return

end subroutine DIAGCT_CPF
