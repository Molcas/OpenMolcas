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

subroutine LOOP70(C,S,ABIJ,AIBJ,AJBI,A,B,F,IPOF,IPOA,IPOB,MYL,NYL,INDA,INDB,INMY,INNY,IFTB,IFTA,FACS,IAB,CPL,CPLA,NVIRA,NVIRC,NVIRB)

use mrci_global, only: NSYM, NVIR
use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: C(*), ABIJ(*), AIBJ(*), AJBI(*), FACS, CPL, CPLA
real(kind=wp), intent(inout) :: S(*)
real(kind=wp), intent(_OUT_) :: A(*), B(*), F(*)
integer(kind=iwp), intent(in) :: IPOF(9), IPOA(9), IPOB(9), MYL, NYL, INDA, INDB, INMY, INNY, IFTB, IFTA
integer(kind=iwp), intent(out) :: IAB, NVIRA, NVIRC, NVIRB
integer(kind=iwp) :: IASYM, IBSYM, ICSYM, IPF
real(kind=wp) :: FACSX

do IASYM=1,NSYM
  IAB = IPOF(IASYM+1)-IPOF(IASYM)
  if (IAB == 0) cycle
  ICSYM = MUL(MYL,IASYM)
  IBSYM = MUL(NYL,ICSYM)
  if ((INDA == INDB) .and. (IBSYM > IASYM)) cycle
  NVIRC = NVIR(ICSYM)
  if (NVIRC == 0) cycle
  NVIRA = NVIR(IASYM)
  NVIRB = NVIR(IBSYM)
  if (ICSYM < IASYM) then
    if (ICSYM < IBSYM) then
      ! CASE 1, IASYM > ICSYM AND IBSYM > ICSYM
      IPF = IPOF(IASYM)+1
      F(1:IAB) = CPL*AIBJ(IPF:IPF+IAB-1)+CPLA*ABIJ(IPF:IPF+IAB-1)
      if (INDA == INDB) call DCOPY_(NVIRA,[Zero],0,F,NVIRA+1)
      call DGEMM_('N','N',NVIRC,NVIRB,NVIRA,FACS,C(INMY+IPOA(IASYM)),NVIRC,F,NVIRA,One,S(INNY+IPOB(IBSYM)),NVIRC)
      if (INDA /= INDB) then
        IPF = IPOF(IBSYM)+1
        F(1:IAB) = CPL*AJBI(IPF:IPF+IAB-1)+CPLA*ABIJ(IPF:IPF+IAB-1)
        call DGEMM_('N','N',NVIRC,NVIRA,NVIRB,FACS,C(INNY+IPOB(IBSYM)),NVIRC,F,NVIRB,One,S(INMY+IPOA(IASYM)),NVIRC)
      end if
    else
      ! CASE 2, IASYM > ICSYM AND ICSYM > OR = IBSYM
      IPF = IPOF(IBSYM)+1
      F(1:IAB) = CPL*AJBI(IPF:IPF+IAB-1)+CPLA*ABIJ(IPF:IPF+IAB-1)

      if (NYL == 1) then
        call DGEMM_('N','T',NVIRB,NVIRC,NVIRA,FACS,F,NVIRB,C(INMY+IPOA(IASYM)),NVIRC,Zero,A,NVIRB)
        if (IFTB == 1) then
          call TRADD(A,S(INNY+IPOB(ICSYM)),NVIRB)
          call SQUARN(C(INNY+IPOB(IBSYM)),A,NVIRB)
        else
          call SIADD(A,S(INNY+IPOB(ICSYM)),NVIRB)
          call SQUAR(C(INNY+IPOB(IBSYM)),A,NVIRB)
        end if
        call DGEMM_('N','N',NVIRC,NVIRA,NVIRB,FACS,A,NVIRC,F,NVIRB,One,S(INMY+IPOA(IASYM)),NVIRC)
      else
        FACSX = FACS
        if (IFTB == 1) FACSX = -FACS
        call DGEMM_('N','T',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,C(INMY+IPOA(IASYM)),NVIRC,One,S(INNY+IPOB(ICSYM)),NVIRB)
        call DGEMM_('T','N',NVIRC,NVIRA,NVIRB,FACSX,C(INNY+IPOB(ICSYM)),NVIRB,F,NVIRB,One,S(INMY+IPOA(IASYM)),NVIRC)
      end if
    end if
  else
    ! UPDATED UNTIL HERE
    if (ICSYM < IBSYM) then
      ! CASE 3, ICSYM > OR = IASYM AND IBSYM > ICSYM
      IPF = IPOF(IASYM)+1
      F(1:IAB) = CPL*AIBJ(IPF:IPF+IAB-1)+CPLA*ABIJ(IPF:IPF+IAB-1)
      if (MYL == 1) then
        if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
        if (IFTA == 1) call SQUARN(C(INMY+IPOA(IASYM)),A,NVIRA)
        call DGEMM_('N','N',NVIRC,NVIRB,NVIRA,FACS,A,NVIRC,F,NVIRA,One,S(INNY+IPOB(IBSYM)),NVIRC)
        call DGEMM_('N','T',NVIRA,NVIRC,NVIRB,FACS,F,NVIRA,C(INNY+IPOB(IBSYM)),NVIRC,Zero,A,NVIRA)
        if (IFTA == 0) call SIADD(A,S(INMY+IPOA(IASYM)),NVIRA)
        if (IFTA == 1) call TRADD(A,S(INMY+IPOA(IASYM)),NVIRA)
      else
        FACSX = FACS
        if (IFTA == 1) FACSX = -FACS
        call DGEMM_('T','N',NVIRC,NVIRB,NVIRA,FACSX,C(INMY+IPOA(ICSYM)),NVIRA,F,NVIRA,One,S(INNY+IPOB(IBSYM)),NVIRC)
        call DGEMM_('N','T',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,C(INNY+IPOB(IBSYM)),NVIRC,One,S(INMY+IPOA(ICSYM)),NVIRA)
      end if
    else
      ! CASE 4, ICSYM > OR = IASYM AND ICSYM > OR = IBSYM
      IPF = IPOF(IBSYM)+1
      F(1:IAB) = CPL*AJBI(IPF:IPF+IAB-1)+CPLA*ABIJ(IPF:IPF+IAB-1)
      if (INDA == INDB) call DCOPY_(NVIRA,[Zero],0,F,NVIRA+1)
      if ((MYL == 1) .and. (NYL == 1)) then

        if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
        if (IFTA == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIRA)
        call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACS,F,NVIRB,A,NVIRA,Zero,B,NVIRB)
        if (IFTB == 0) call SIADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
        if (IFTB == 1) call TRADD(B,S(INNY+IPOB(ICSYM)),NVIRB)

      else if ((MYL == 1) .and. (NYL /= 1)) then
        if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
        if (IFTA == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIRA)
        FACSX = FACS
        if (IFTB == 1) FACSX = -FACS
        call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,A,NVIRA,One,S(INNY+IPOB(ICSYM)),NVIRB)

      else if ((MYL /= 1) .and. (NYL == 1)) then

        FACSX = FACS
        if (IFTA == 1) FACSX = -FACS
        call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,C(INMY+IPOA(ICSYM)),NVIRA,Zero,B,NVIRB)
        if (IFTB == 0) call SIADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
        if (IFTB == 1) call TRADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
      else if ((MYL /= 1) .and. (NYL /= 1)) then
        FACSX = FACS
        if (IFTA+IFTB == 1) FACSX = -FACS
        call DGEMM_('N','N',NVIRB,NVIRC,NVIRA,FACSX,F,NVIRB,C(INMY+IPOA(ICSYM)),NVIRA,One,S(INNY+IPOB(ICSYM)),NVIRB)
      end if
      if (INDA /= INDB) then
        IPF = IPOF(IASYM)+1
        F(1:IAB) = CPL*AIBJ(IPF:IPF+IAB-1)+CPLA*ABIJ(IPF:IPF+IAB-1)

        if ((NYL == 1) .and. (MYL == 1)) then

          if (IFTB == 0) call SQUAR(C(INNY+IPOB(IBSYM)),A,NVIRB)
          if (IFTB == 1) call SQUARM(C(INNY+IPOB(IBSYM)),A,NVIRB)
          call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACS,F,NVIRA,A,NVIRB,Zero,B,NVIRA)
          if (IFTA == 0) call SIADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
          if (IFTA == 1) call TRADD(B,S(INMY+IPOA(ICSYM)),NVIRA)

        else if ((NYL == 1) .and. (MYL /= 1)) then

          if (IFTB == 0) call SQUAR(C(INNY+IPOB(ICSYM)),A,NVIRB)
          if (IFTB == 1) call SQUARM(C(INNY+IPOB(ICSYM)),A,NVIRB)
          FACSX = FACS
          if (IFTA == 1) FACSX = -FACS
          call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,A,NVIRB,One,S(INMY+IPOA(ICSYM)),NVIRA)

        else if ((NYL /= 1) .and. (MYL == 1)) then

          FACSX = FACS
          if (IFTB == 1) FACSX = -FACS
          call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,C(INNY+IPOB(ICSYM)),NVIRB,Zero,B,NVIRA)
          if (IFTA == 0) call SIADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
          if (IFTA == 1) call TRADD(B,S(INMY+IPOA(ICSYM)),NVIRA)

        else if ((NYL /= 1) .and. (MYL /= 1)) then

          FACSX = FACS
          if (IFTA+IFTB == 1) FACSX = -FACS
          call DGEMM_('N','N',NVIRA,NVIRC,NVIRB,FACSX,F,NVIRA,C(INNY+IPOB(ICSYM)),NVIRB,One,S(INMY+IPOA(ICSYM)),NVIRA)

        end if
      end if
    end if
  end if

end do

return

end subroutine LOOP70
