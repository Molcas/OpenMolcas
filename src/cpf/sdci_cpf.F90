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

subroutine SDCI_CPF(MEMORY)

use cpf_global, only: ICASE, ICONV, ICPF, IDENS, ILIM, INCPF, INDX, IPRINT, IR1, IRC, IREF0, IREST, IROW, ISDCI, ISMAX, ITER, &
                      ITPUL, JMAX, JSC, JSY, KBUFF1, LBUF, LIC, LN, MAXIT, MX1, MX2, NORBT, NREF, NTMAX, NVMAX
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6, RtoI

implicit none
integer(kind=iwp), intent(in) :: MEMORY
integer(kind=iwp), allocatable :: IBMN(:), ICASE_(:)
real(kind=wp), allocatable :: A(:), ABIJ(:), AC1(:), AC2(:), AIBJ(:), AJBI(:), AP(:), B(:), BIJ(:), BMN(:), BST(:), BUFAC(:), &
                              BUFIJ(:), BUFIN(:), C(:), CN(:), DBK(:), ENP(:), EPB(:), EPP(:), F(:), FC(:), FIJKL(:), FK(:), &
                              FSEC(:), S(:), TEMP(:), TEMP2(:), THET(:), TPQ(:), W(:)

LIC = MEMORY
IPRINT = 5
IDENS = 0
! INPUT, SORTING AND DIAGONAL ELEMENTS
call READIN_CPF()
call DIAGCT_CPF()
ITER = 1
if (IREST == 1) ITER = ITER+1
ITPUL = 1
call mma_allocate(C,JSC(ILIM),label='C')
if (IREST == 0) call START_CPF(C,JSC(4),IREF0)
if (IREST == 1) call RESTART_CPFMCPF(C,JSC(4))
call mma_allocate(THET,IRC(ILIM)**2,label='THET')
if ((ICPF == 0) .and. (INCPF == 0) .and. (ISDCI == 0)) then
  call THETSET(ICASE,THET,IRC(4))
  call mma_allocate(W,JSC(ILIM),label='W')
else
  call mma_allocate(W,0,label='W')
end if
call mma_allocate(FC,IROW(NORBT+1),label='FC')
call mma_allocate(BUFIN,LBUF+(LBUF+2+(RtoI-1))/RtoI,label='BUFIN') ! LBUF reals + LBUF+2 integers
call mma_allocate(A,MX2,label='A')
call mma_allocate(B,MX2,label='B')
call mma_allocate(FK,max(NVMAX,MX1),label='FK')
call mma_allocate(DBK,NVMAX,label='DBK')
call mma_allocate(TEMP,IRC(ILIM),label='TEMP')
call mma_allocate(S,JSC(ILIM),label='S')
call mma_allocate(TPQ,IRC(ILIM),label='TPQ')
call mma_allocate(ENP,IRC(ILIM),label='ENP')
call mma_allocate(EPP,IRC(ILIM),label='EPP')
call mma_allocate(BST,(MAXIT+1)**2,label='BST')
call mma_allocate(ABIJ,MX1,label='ABIJ')
call mma_allocate(AIBJ,MX1,label='AIBJ')
call mma_allocate(AJBI,MX1,label='AJBI')
call mma_allocate(F,MX1,label='F')
call mma_allocate(FSEC,2*MX1,label='FSEC') ! it was MX1, but that's not enough
call mma_allocate(FIJKL,(IROW(LN+1)*(IROW(LN+1)+1))/2,label='FIJKL')
call mma_allocate(BUFIJ,KBUFF1+2,label='BUFIJ')
call mma_allocate(BMN,JMAX,label='BMN')
call mma_allocate(IBMN,JMAX,label='IBMN')
call mma_allocate(AC1,ISMAX,label='AC1')
call mma_allocate(AC2,ISMAX,label='AC2')
call mma_allocate(BUFAC,KBUFF1,label='BUFAC')
call mma_allocate(EPB,IRC(ILIM),label='EPB')
call mma_allocate(AP,IRC(ILIM),label='AP')
call mma_allocate(BIJ,(MAXIT+1)**2,label='BIJ')
call mma_allocate(CN,MAXIT+1,label='CN')
call mma_allocate(TEMP2,NTMAX,label='TEMP2')
do
  call NPSET(JSY,INDX,C,TPQ,ENP,TEMP,S,W,EPP,ICASE)
  call TWOCT(C,S,W,THET,ENP,EPP,ABIJ,AIBJ,AJBI,BUFIN,A,B,F,FSEC,FIJKL,BUFIJ,BMN,IBMN,AC1,AC2,BUFAC)
  call ONECT(C,S,W,THET,ENP,EPP,FC,BUFIN,A,B,FK,DBK)
  call CPFCTL(C,S,W,TPQ,ENP,EPP,BST,EPB,AP,BIJ,CN,TEMP2)
  ITER = ITER+1
  ITPUL = ITPUL+1
  if ((ITER > MAXIT) .or. (ICONV == 1)) exit
end do
call mma_deallocate(BST)
call mma_deallocate(ABIJ)
call mma_deallocate(AIBJ)
call mma_deallocate(AJBI)
call mma_deallocate(F)
call mma_deallocate(FSEC)
call mma_deallocate(FIJKL)
call mma_deallocate(BUFIJ)
call mma_deallocate(BMN)
call mma_deallocate(IBMN)
call mma_deallocate(AC1)
call mma_deallocate(AC2)
call mma_deallocate(BUFAC)
call mma_deallocate(EPB)
call mma_deallocate(AP)
call mma_deallocate(BIJ)
call mma_deallocate(CN)
call mma_deallocate(TEMP2)
! FIRST ORDER DENSITY MATRIX
IDENS = 1

call mma_allocate(ICASE_,IR1,label='ICASE')
call DENSCT_CPF(C,S,W,THET,TPQ,ENP,EPP,ICASE_,FC,BUFIN,A,B,FK,DBK,TEMP)
call mma_deallocate(C)
call mma_deallocate(S)
call mma_deallocate(W)
call mma_deallocate(THET)
call mma_deallocate(TPQ)
call mma_deallocate(ENP)
call mma_deallocate(EPP)
call mma_deallocate(ICASE_)
call mma_deallocate(FC)
call mma_deallocate(BUFIN)
call mma_deallocate(A)
call mma_deallocate(B)
call mma_deallocate(FK)
call mma_deallocate(DBK)
call mma_deallocate(TEMP)

if (NREF > 1) then
  write(u6,*) ' This is a single reference program, but more than'
  write(u6,*) ' one reference state has been specified in the'
  write(u6,*) ' GUGA program. Change input to GUGA and run again.'
end if

return

end subroutine SDCI_CPF
