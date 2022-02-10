!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine LOOP70(INTSYM,INDX,C,S,ABIJ,AIBJ,AJBI,BUF,IBUF,A,B,F,FSEC,IPOF,IPOA,IPOB,MYL,NYL,INDA,INDB,INMY,INNY,IFTB,IFTA,FACS, &
                  IAB,CPL,CPLA,NVIRA,NVIRC,NVIRB)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
dimension INTSYM(*), INDX(*), C(*), S(*), ABIJ(*), AIBJ(*), AJBI(*), BUF(NBITM3), IBUF(NBITM3+2), A(*), B(*), F(*), FSEC(*)
dimension IPOF(9), IPOA(9), IPOB(9)

do IASYM=1,NSYM
  IAB = IPOF(IASYM+1)-IPOF(IASYM)
  if (IAB == 0) GO TO 70
  ICSYM = MUL(MYL,IASYM)
  IBSYM = MUL(NYL,ICSYM)
  if ((INDA == INDB) .and. (IBSYM > IASYM)) GO TO 70
  NVIRC = NVIR(ICSYM)
  if (NVIRC == 0) GO TO 70
  NVIRA = NVIR(IASYM)
  NVIRB = NVIR(IBSYM)
  if (ICSYM >= IASYM) GO TO 31
  if (ICSYM >= IBSYM) GO TO 32
  ! CASE 1, IASYM > ICSYM AND IBSYM > ICSYM
  IPF = IPOF(IASYM)+1
  call DYAX(IAB,CPL,AIBJ(IPF),1,F,1)
  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
  if (INDA == INDB) call SETZZ(F,NVIRA)
  !call FZERO(A,NBC)
  !call FMMM(C(INMY+IPOA(IASYM)),F,A,NVIRC,NVIRB,NVIRA)
  !call DAXPY_(NBC,FACS,A,1,S(INNY+IPOB(IBSYM)),1)
  call DGEMM_('N','N',NVIRC,NVIRB,NVIRA,FACS,C(INMY+IPOA(IASYM)),NVIRC,F,NVIRA,1.0d00,S(INNY+IPOB(IBSYM)),NVIRC)
  if (INDA == INDB) GO TO 70
  IPF = IPOF(IBSYM)+1
  call FZERO(F,IAB)
  call DYAX(IAB,CPL,AJBI(IPF),1,F,1)
  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
  call DGEMM_('N','N',NVIRC,NVIRA,NVIRB,FACS,C(INNY+IPOB(IBSYM)),NVIRC,F,NVIRB,1.0d00,S(INMY+IPOA(IASYM)),NVIRC)
  GO TO 70
  ! CASE 2, IASYM > ICSYM AND ICSYM > OR = IBSYM
32 IPF = IPOF(IBSYM)+1
  call DYAX(IAB,CPL,AJBI(IPF),1,F,1)
  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)

  if (NYL == 1) then
    call DGEMM_('N','T',NVIRB,NVIRC,NVIRA,FACS,F,NVIRB,C(INMY+IPOA(IASYM)),NVIRC,0.d0,A,NVIRB)
    if (IFTB == 1) then
      call TRADD(A,S(INNY+IPOB(ICSYM)),NVIRB)
      call SQUARN(C(INNY+IPOB(IBSYM)),A,NVIRB)
    else
      call SIADD(A,S(INNY+IPOB(ICSYM)),NVIRB)
      call SQUAR(C(INNY+IPOB(IBSYM)),A,NVIRB)
    end if
    call DGEMM_('N','N',NVIRC,NVIRA,NVIRB,FACS,A,NVIRC,F,NVIRB,1.d0,S(INMY+IPOA(IASYM)),NVIRC)
  else
    FACSX = FACS
    if (IFTB == 1) FACSX = -FACS
    call DGEMM_('N','T',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,C(INMY+IPOA(IASYM)),NVIRC,1.d0,S(INNY+IPOB(ICSYM)),NVIRB)
    call DGEMM_('T','N',NVIRC,NVIRA,NVIRB,FACSX,C(INNY+IPOB(ICSYM)),NVIRB,F,NVIRB,1.d0,S(INMY+IPOA(IASYM)),NVIRC)
  end if
  GO TO 70
  ! UPDATED UNTIL HERE
31 if (ICSYM >= IBSYM) GO TO 33
  ! CASE 3, ICSYM > OR = IASYM AND IBSYM > ICSYM
  IPF = IPOF(IASYM)+1
  call DYAX(IAB,CPL,AIBJ(IPF),1,F,1)
  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
  if (MYL == 1) then
    if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
    if (IFTA == 1) call SQUARN(C(INMY+IPOA(IASYM)),A,NVIRA)
    call DGEMM_('N','N',NVIRC,NVIRB,NVIRA,FACS,A,NVIRC,F,NVIRA,1.d0,S(INNY+IPOB(IBSYM)),NVIRC)
    call DGEMM_('N','T',NVIRA,NVIRC,NVIRB,FACS,F,NVIRA,C(INNY+IPOB(IBSYM)),NVIRC,0.d0,A,NVIRA)
    if (IFTA == 0) call SIADD(A,S(INMY+IPOA(IASYM)),NVIRA)
    if (IFTA == 1) call TRADD(A,S(INMY+IPOA(IASYM)),NVIRA)
  else
    FACSX = FACS
    if (IFTA == 1) FACSX = -FACS
    call DGEMM_('T','N',NVIRC,NVIRB,NVIRA,FACSX,C(INMY+IPOA(ICSYM)),NVIRA,F,NVIRA,1.d0,S(INNY+IPOB(IBSYM)),NVIRC)
    call DGEMM_('N','T',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,C(INNY+IPOB(IBSYM)),NVIRC,1.d0,S(INMY+IPOA(ICSYM)),NVIRA)
  end if
  GO TO 70
  ! CASE 4, ICSYM > OR = IASYM AND ICSYM > OR = IBSYM
33 IPF = IPOF(IBSYM)+1
  call DYAX(IAB,CPL,AJBI(IPF),1,F,1)
  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
  if (INDA == INDB) call SETZZ(F,NVIRA)
  if ((MYL == 1) .and. (NYL == 1)) then

    if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
    if (IFTA == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIRA)
    call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACS,F,NVIRB,A,NVIRA,0.d0,B,NVIRB)
    if (IFTB == 0) call SIADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
    if (IFTB == 1) call TRADD(B,S(INNY+IPOB(ICSYM)),NVIRB)

  else if ((MYL == 1) .and. (NYL /= 1)) then
    if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
    if (IFTA == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIRA)
    FACSX = FACS
    if (IFTB == 1) FACSX = -FACS
    call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,A,NVIRA,1.d0,S(INNY+IPOB(ICSYM)),NVIRB)

  else if ((MYL /= 1) .and. (NYL == 1)) then

    FACSX = FACS
    if (IFTA == 1) FACSX = -FACS
    call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,C(INMY+IPOA(ICSYM)),NVIRA,0.d0,B,NVIRB)
    if (IFTB == 0) call SIADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
    if (IFTB == 1) call TRADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
  else if ((MYL /= 1) .and. (NYL /= 1)) then
    FACSX = FACS
    if (IFTA+IFTB == 1) FACSX = -FACS
    call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,C(INMY+IPOA(ICSYM)),NVIRA,1.d0,S(INNY+IPOB(ICSYM)),NVIRB)
  end if
  if (INDA == INDB) GO TO 70
  IPF = IPOF(IASYM)+1
  call DYAX(IAB,CPL,AIBJ(IPF),1,F,1)
  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)

  if ((NYL == 1) .and. (MYL == 1)) then

    if (IFTB == 0) call SQUAR(C(INNY+IPOB(IBSYM)),A,NVIRB)
    if (IFTB == 1) call SQUARM(C(INNY+IPOB(IBSYM)),A,NVIRB)
    call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACS,F,NVIRA,A,NVIRB,0.d0,B,NVIRA)
    if (IFTA == 0) call SIADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
    if (IFTA == 1) call TRADD(B,S(INMY+IPOA(ICSYM)),NVIRA)

  else if ((NYL == 1) .and. (MYL /= 1)) then
    if (IFTB == 0) call SQUAR(C(INNY+IPOB(ICSYM)),A,NVIRB)
    if (IFTB == 1) call SQUARM(C(INNY+IPOB(ICSYM)),A,NVIRB)
    FACSX = FACS
    if (IFTA == 1) FACSX = -FACS
    call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,A,NVIRB,1.d0,S(INMY+IPOA(ICSYM)),NVIRA)

  else if ((NYL /= 1) .and. (MYL == 1)) then

    FACSX = FACS
    if (IFTB == 1) FACSX = -FACS
    call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,C(INNY+IPOB(ICSYM)),NVIRB,0.d0,B,NVIRA)
    if (IFTA == 0) call SIADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
    if (IFTA == 1) call TRADD(B,S(INMY+IPOA(ICSYM)),NVIRA)

  else if ((NYL /= 1) .and. (MYL /= 1)) then

    FACSX = FACS
    if (IFTA+IFTB == 1) FACSX = -FACS
    call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,C(INNY+IPOB(ICSYM)),NVIRB,1.d0,S(INMY+IPOA(ICSYM)),NVIRA)
  end if

70 continue
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(INTSYM)
  call Unused_integer_array(INDX)
  call Unused_integer_array(IBUF)
  call Unused_real_array(FSEC)
  call Unused_real_array(BUF)
end if

end subroutine LOOP70
