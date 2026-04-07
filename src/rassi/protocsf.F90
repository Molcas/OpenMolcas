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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************

subroutine PROTOCSF(NPEL,MLTPL,NPCSFSZ,IPCSFCP)
! Return a table with all possible CSF's with NPEL electrons
! coupled to give a total spin multiplicity MLTPL.
! The order of the resulting CSF's is consistent with the
! index function,
!  Index=1+Sum(j) NGENE(j-1,2*S_j+2)
! where the sum is over only the up-coupled orbitals j,
! S_j is the accumulated spin, summed over orbitals <= j,
! and NGENE(N,2*S+1) is in general the number of genealogical
! couplings of N electrons to obtain spin S.

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NPEL, MLTPL, NPCSFSZ, IPCSFCP(NPEL,NPCSFSZ)
integer(kind=iwp) :: ISP2, N, ND, NPCSF, NPELD, NPELU, NU
integer(kind=iwp), parameter :: DWNCPL = 0, UPCPL = 1
integer(kind=iwp), external :: NGENE

if (NPEL == 0) return
ISP2 = MLTPL-1
if (ISP2 < 0) return
if (ISP2 > NPEL) return
NPELU = (NPEL+ISP2)/2
NPELD = (NPEL-ISP2)/2
if (NPELU < NPELD) return
if (NPELU+NPELD /= NPEL) return
IPCSFCP(1:NPELU,1) = UPCPL
if (NPELU == NPEL) return
IPCSFCP(NPELU+1:NPEL,1) = DWNCPL

NPCSF = NGENE(NPEL,MLTPL)
if (NPCSF > NPCSFSZ) then
  write(u6,*) ' Too small space allocated in PROTOCSF. Input:'
  write(u6,'(1x,a,3i6)') ' NPEL,MLTPL,NPCSFSZ:',NPEL,MLTPL,NPCSFSZ
  write(u6,'(1x,a,i12)') ' Required NPCSFSZ is',NPCSF
  call ABEND()
end if

NPCSF = 1
if (NPEL <= 2) return
outer: do
  N = 0
  NU = 0
  do
    N = N+1
    if (N > NPEL) exit outer
    if (IPCSFCP(N,NPCSF) == UPCPL) NU = NU+1
    ND = N-NU
    if ((IPCSFCP(N,NPCSF) /= UPCPL) .and. (NU /= ND)) exit
  end do
  IPCSFCP(1:NU-1,NPCSF+1) = UPCPL
  IPCSFCP(NU:N-1,NPCSF+1) = DWNCPL
  IPCSFCP(N,NPCSF+1) = UPCPL
  IPCSFCP(N+1:NPEL,NPCSF+1) = IPCSFCP(N+1:NPEL,NPCSF)
  NPCSF = NPCSF+1
end do outer

end subroutine PROTOCSF
