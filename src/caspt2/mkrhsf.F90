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
!--------------------------------------------
! 1998  PER-AAKE MALMQUIST
! DEPARTMENT OF THEORETICAL CHEMISTRY
! UNIVERSITY OF LUND
! SWEDEN
!--------------------------------------------

subroutine MKRHSF(IVEC,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for cases 8 and 9 (BVAT).

use Symmetry_Info, only: Mul
use SUPERINDEX, only: KAGEB, KAGTB, KTGEU, KTGTU
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array, GA_Arrays
use general_data, only: NASH
use caspt2_module, only: NAES, NAGEBES, NAGTBES, NASUP, NINDEP, NISH, NISUP, NORB, NSES, NSSH, NSYM, NTGEUES, NTGTUES
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, nERI1, nERI2, nSCR
real(kind=wp), intent(inout) :: ERI1(nERI1), ERI2(nERI2), SCR(nSCR)
integer(kind=iwp) :: IA, IAABS, IATOT, IB, IBABS, IBTOT, IBUF, ICASE, ISYM, ISYMA, ISYMB, ISYMT, ISYMU, IT, ITABS, ITTOT, IU, &
                     IUABS, IUTOT, IWAM, IWAP, IWIM, IWIP, IWM, JWP, LWM, LWP, NASM, NASP, NINM, NINP, NISM, NISP, NVM, NVP
real(kind=wp) :: A, B
real(kind=wp), parameter :: SQI2 = sqrt(Half)

do ISYM=1,NSYM
  NINP = NINDEP(ISYM,8)
  NINM = NINDEP(ISYM,9)
  if (NINP+NINM == 0) cycle
  NASP = NASUP(ISYM,8)
  NISP = NISUP(ISYM,8)
  NASM = NASUP(ISYM,9)
  NISM = NISUP(ISYM,9)
  NVP = NASP*NISP
  if (NVP == 0) cycle
  NVM = NASM*NISM
  LWP = Allocate_GA_Array(NVP,'WFP')
  if (NVM > 0) LWM = Allocate_GA_Array(NVM,'WFM')
  !  Let W(t,u,ab)=(aubt)
  !   WP(tu,ab)=(W(t,u,ab)+W(u,t,ab))*(1-Kron(t,u)/2) /2
  ! With new normalisation, replace /2 with /(2*SQRT(1+Kron(ab))
  !   WM(tu,ab)=(W(t,u,ab)-W(u,t,ab))*(1-Kron(t,u)/2) /2
  do ISYMA=1,NSYM
    ISYMB = Mul(ISYMA,ISYM)
    if (ISYMA < ISYMB) cycle
    do ISYMT=1,NSYM
      ISYMU = Mul(ISYMT,ISYM)
      if (ISYMT < ISYMU) cycle
      do IT=1,NASH(ISYMT)
        ITABS = IT+NAES(ISYMT)
        ITTOT = IT+NISH(ISYMT)
        do IU=1,NASH(ISYMU)
          IUABS = IU+NAES(ISYMU)
          IUTOT = IU+NISH(ISYMU)
          if (ITABS < IUABS) exit
          call EXCH(ISYMA,ISYMU,ISYMB,ISYMT,IUTOT,ITTOT,ERI1,SCR)
          call EXCH(ISYMA,ISYMT,ISYMB,ISYMU,ITTOT,IUTOT,ERI2,SCR)
          do IA=1,NSSH(ISYMA)
            IAABS = IA+NSES(ISYMA)
            IATOT = IA+NISH(ISYMA)+NASH(ISYMA)
            do IB=1,NSSH(ISYMB)
              IBABS = IB+NSES(ISYMB)
              IBTOT = IB+NISH(ISYMB)+NASH(ISYMB)
              if (IAABS < IBABS) exit
              IBUF = IATOT+NORB(ISYMA)*(IBTOT-1)
              A = Half*(ERI1(IBUF)+ERI2(IBUF))
              if (ITABS == IUABS) A = Half*A
              IWAP = KTGEU(ITABS,IUABS)-NTGEUES(ISYM)
              IWIP = KAGEB(IAABS,IBABS)-NAGEBES(ISYM)
              JWP = IWAP+NASP*(IWIP-1)
              if (IAABS /= IBABS) then
                GA_Arrays(LWP)%A(JWP) = A
                if (ITABS /= IUABS) then
                  B = Half*(ERI1(IBUF)-ERI2(IBUF))
                  IWAM = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
                  IWIM = KAGTB(IAABS,IBABS)-NAGTBES(ISYM)
                  IWM = IWAM+NASM*(IWIM-1)
                  GA_Arrays(LWM)%A(IWM) = B
                end if
              else
                GA_Arrays(LWP)%A(JWP) = SQI2*A
              end if
            end do
          end do
        end do
      end do
    end do
  end do
  ! Put WP on disk
  ICASE = 8
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
  call Deallocate_GA_Array(LWP)
  if (NINM > 0) then
    ! Put WM on disk
    ICASE = 9
    call MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
  end if
  if (NVM > 0) call Deallocate_GA_Array(LWM)
end do

end subroutine MKRHSF
