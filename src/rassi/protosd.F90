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

subroutine ProtoSD(NPELA,NPELB,NPSDSZ,IPSDMS)
! Given NPELA and NPELB, returns a table of all possible
! Slater determinants with NPELA alpha electrons and NPELB
! beta electrons in NPELA+NPELB orbitals.
! The ordering of the resulting SD is consistent with the
! index function   Index=1+Sum(j) noverm(k-1,j), where the sum
! is over only the alpha spins, enumerated j=1..NPELA,
! while k is the orbital with the j-th alpha electron.

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NPELA, NPELB, NPSDSZ, IPSDMS(NPELA+NPELB,NPSDSZ)
integer(kind=iwp) :: ITMP(50), J, K, L, NDET, NPORB
integer(kind=iwp), parameter :: ASPIN = 1, BSPIN = 0
integer(kind=iwp), external :: NOVERM

if (NPELA < 0) goto 999
if (NPELB < 0) goto 999
NPORB = NPELA+NPELB
if (NPORB == 0) return
do K=1,NPELA
  ITMP(K) = K
  IPSDMS(K,1) = ASPIN
end do
if (NPELA == NPORB) return
do K=NPELA+1,NPORB
  IPSDMS(K,1) = BSPIN
end do
if (NPELA == 0) return

NDET = NOVERM(NPORB,NPELA)
if (NDET > NPSDSZ) goto 998
!PAM write(u6,*) 'PROTOSD. Nr of orbitals:',NPORB
!PAM write(u6,*) '         Nr of determin:',NDET
!PAM write(u6,*) '         Need array siz:',NPORB*NDET

ITMP(NPELA+1) = NPORB+1
NDET = 1

10 continue
K = 0
20 continue
K = K+1
if (K > NPELA) goto 30
if (ITMP(K+1) == 1+ITMP(K)) goto 20
ITMP(K) = ITMP(K)+1
do L=1,K-1
  ITMP(L) = L
end do
NDET = NDET+1
if (NDET > NPSDSZ) goto 997
do J=1,NPORB
  IPSDMS(J,NDET) = BSPIN
end do
do L=1,NPELA
  IPSDMS(ITMP(L),NDET) = ASPIN
end do
goto 10

30 continue
return
997 continue
write(u6,*) " Serious error in PROTOSD. Too many SD's are produced."
call ABEND()
998 continue
write(u6,*) ' Too small space allocated in PROTOSD. Input:'
write(u6,'(1x,a,3i6)') ' NPELA,NPELB,NPSDSZ:',NPELA,NPELB,NPSDSZ
write(u6,'(1x,a,i12)') ' Required NPSDSZ is',NDET
call ABEND()
999 continue
write(u6,*) ' Invalid input to ProtoSD.'
write(u6,*) '  NPELA,NPELB:',NPELA,NPELB
call ABEND()

end subroutine ProtoSD
