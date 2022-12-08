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
!***********************************************************************
!***********************************************************************
!                                                                      *
! PER SIEGBAHN                                                         *
! DEPARTMENT OF THEORETICAL PHYSICS                                    *
! UNIVERSITY OF STOCKHOLM                                              *
! SWEDEN                                                               *
!                                                                      *
! UPDATED FOR MOLCAS-4 BY P-A MALMQVIST & NIGEL MORIARTY 1996          *
!***********************************************************************

subroutine GUGA(IRETURN)

use guga_global, only: free_all, IADD10, IPRINT, Lu_10, Lu_11, MXVERT, NBUF
use guga_util_global, only: IAD10
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Four
use Definitions, only: wp, iwp, u6, RtoI

implicit none
integer(kind=iwp), intent(out) :: IRETURN
integer(kind=iwp) :: ISPAC, IST, KB, KBUF, KBUF2, LSTO, LW1, MCOP, NBINS, NCOR, NCORX, NTPB
real(kind=wp) :: A, B, C
integer(kind=iwp), allocatable :: L0(:), L1(:), L2(:), L3(:)
integer(kind=iwp), external :: isFreeUnit

! Prologue

! Find max possible allocatable
call mma_maxINT(NCOR)
! Grab almost all of it, but leave a little to be safe:
NCOR = NCOR-100000
!call SetQue('Trace=On')
call SETTIM()
!call XUFLOW()
!call ERRSET(208,256,-1,1,1,208)
!call ERRSET(151,256,-1,1,1,151)
call JTIME(IST)

! Print program header

! Initialize files

Lu_11 = isFreeUnit(11)
call DANAME_wa(Lu_11,'TEMP01')
Lu_10 = isFreeUnit(10)
call DANAME(Lu_10,'CIGUGA')
IAD10(:) = 0
IADD10 = 0
call iDAFILE(Lu_10,1,IAD10,9,IADD10)

! Read input

NBUF = 600
!PAM96: Use variable MCOP, size of buffers:
MCOP = NBUF*RtoI+NBUF+1
call mma_allocate(L0,4*MXVERT,label='L0')
call mma_allocate(L1,4*MXVERT,label='L1')
call mma_allocate(L2,4*MXVERT,label='L2')
call mma_allocate(L3,4*MXVERT,label='L3')
call INPUT_GUGA(L0,L1,L2,L3,ISPAC)

! Main body

! SORT ALLOCATION , ISPAC WORDS TO BE SORTED IN NCOR CORE SPACE
! TWO BUFFERS OF LENGTH NBINS NEEDED
! NCOR IS IN UNITS OF FLOATING-POINT WORDS, e.g. REAL*8
! NBINS=ISPAC/(NCORX-2*NBINS)+1
NCORX = NCOR-(MCOP+1)
A = Two
B = -NCORX-2
C = ISPAC+NCORX
NBINS = int((-B-sqrt(B*B-Four*A*C))/(Two*A))
! NUMBER OF WORDS IN EACH BIN
NTPB = (ISPAC-1)/NBINS+1
! SPACE IN CORE FOR EACH BIN
KB = RtoI*(NCORX-2*NBINS)/NBINS
KBUF = (KB-1)/(RtoI+1)
KBUF = (KBUF/2)*2
if (KBUF > 600) KBUF = 600
KBUF2 = KBUF*RtoI+KBUF+2
if (IPRINT >= 2) write(u6,10) KBUF,NBINS,NTPB,NCOR,ISPAC
! STORAGE FOR NBINS BINS EACH OF SIZE KBUF2 IN AIAI
LSTO = NBINS*KBUF2
! ALSO SPACE FOR NTPB WORDS IN EMPTY
if (NTPB > LSTO) LSTO = NTPB
LW1 = LSTO+1
call CI_SELECT(L0,L1,L2,L3,KBUF,NTPB,NBINS,LW1)
call mma_deallocate(L0)
call mma_deallocate(L1)
call mma_deallocate(L2)
call mma_deallocate(L3)
IADD10 = 0
call iDAFILE(Lu_10,1,IAD10,9,IADD10)

! Epilogue, end
call free_all()

!                                                                      *
!***********************************************************************
!                                                                      *
! Close open dafiles

call DaClos(Lu_10)
call DaClos(Lu_11)

!                                                                      *
!***********************************************************************
!                                                                      *
ireturn = 0

return

10 format(/6X,'SORTING INFORMATION',/6X,'KBUF=',I7,/6X,'NBINS=',I6,/6X,'NTPB=',I7,/6X,'NCOR=',I7,/6X,'ISPAC=',I6)

end subroutine GUGA
