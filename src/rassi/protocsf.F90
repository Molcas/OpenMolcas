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

use Definitions, only: u6

dimension IPCSFCP(NPEL,NPCSFSZ)
integer UPCPL, DWNCPL
parameter(UPCPL=1,DWNCPL=0)

if (NPEL == 0) return
ISP2 = MLTPL-1
if (ISP2 < 0) return
if (ISP2 > NPEL) return
NPELU = (NPEL+ISP2)/2
NPELD = (NPEL-ISP2)/2
if (NPELU < NPELD) return
if (NPELU+NPELD /= NPEL) return
do N=1,NPELU
  IPCSFCP(N,1) = UPCPL
end do
if (NPELU == NPEL) return
do N=NPELU+1,NPEL
  IPCSFCP(N,1) = DWNCPL
end do

NPCSF = NGENE(NPEL,MLTPL)
if (NPCSF > NPCSFSZ) goto 998

NPCSF = 1
if (NPEL <= 2) return
10 continue
N = 0
NU = 0
20 continue
N = N+1
if (N > NPEL) goto 30
if (IPCSFCP(N,NPCSF) == UPCPL) NU = NU+1
ND = N-NU
if ((IPCSFCP(N,NPCSF) == UPCPL) .or. (NU == ND)) goto 20
do L=1,NU-1
  IPCSFCP(L,NPCSF+1) = UPCPL
end do
do L=NU,N-1
  IPCSFCP(L,NPCSF+1) = DWNCPL
end do
IPCSFCP(N,NPCSF+1) = UPCPL
do L=N+1,NPEL
  IPCSFCP(L,NPCSF+1) = IPCSFCP(L,NPCSF)
end do
NPCSF = NPCSF+1
goto 10

30 continue
return
998 continue
write(u6,*) ' Too small space allocated in PROTOCSF. Input:'
write(u6,'(1x,a,3i6)') ' NPEL,MLTPL,NPCSFSZ:',NPEL,MLTPL,NPCSFSZ
write(u6,'(1x,a,i12)') ' Required NPCSFSZ is',NPCSF
call ABEND()

end subroutine PROTOCSF
