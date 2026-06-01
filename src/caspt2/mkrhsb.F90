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

subroutine MKRHSB(IVEC,ERI,nERI,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV for cases 2 and 3 (VJTI).

use Symmetry_Info, only: Mul
use SUPERINDEX, only: KIGEJ, KIGTJ, KTGEU, KTGTU
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array, GA_Arrays
use caspt2_module, only: NAES, NASH, NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NINDEP, NISH, NORB, NSYM, NTGEU, NTGEUES, NTGTU, NTGTUES
use Constants, only: Two, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, nERI, nSCR
real(kind=wp), intent(inout) :: ERI(nERI), SCR(nSCR)
integer(kind=iwp) :: IBUF, ICASE, II, IIABS, IIJM, IIJP, IJ, IJABS, ISYM, ISYMI, ISYMJ, ISYMT, ISYMU, IT, ITABS, ITTOT, ITUM, &
                     ITUP, IU, IUABS, IUTOT, IWM, JWP, LWM, LWP, NASM, NASP, NINM, NINP, NISM, NISP, NVM, NVP
real(kind=wp) :: Val
real(kind=wp), parameter :: SQ2 = sqrt(Two)

! VJTI CASE:
do ISYM=1,NSYM
  NINP = NINDEP(ISYM,2)
  NINM = NINDEP(ISYM,3)
  if (NINP+NINM == 0) cycle
  NASP = NTGEU(ISYM)
  NISP = NIGEJ(ISYM)
  NVP = NASP*NISP
  if (NVP == 0) cycle
  NASM = NTGTU(ISYM)
  NISM = NIGTJ(ISYM)
  NVM = NASM*NISM
  ! Allocate WP,WM
  LWP = Allocate_GA_Array(NVP,'WBP')
  LWM = Allocate_GA_Array(NVM,'WBM')
  ! Let  W(tu,i,j)=(it,ju):
  !   WP(tu,ij)=(W(tu,i,j)+W(tu,j,i))*(1-Kron(t,u)/2) /2
  ! With new normalisation, replace /2 with /(2*SQRT(1+Kron(ij))
  !   WM(tu,ij)=(W(tu,i,j)-W(tu,j,i))*(1-Kron(t,u)/2) /2
  do ISYMT=1,NSYM
    ISYMU = Mul(ISYMT,ISYM)
    if (ISYMT < ISYMU) cycle
    if (NASH(ISYMT)*NASH(ISYMU) == 0) cycle
    do ISYMI=1,NSYM
      ISYMJ = Mul(ISYMI,ISYM)
      if (NISH(ISYMI)*NISH(ISYMJ) == 0) cycle
      do IT=1,NASH(ISYMT)
        ITABS = IT+NAES(ISYMT)
        ITTOT = IT+NISH(ISYMT)
        do IU=1,NASH(ISYMU)
          IUABS = IU+NAES(ISYMU)
          IUTOT = IU+NISH(ISYMU)
          if (ITABS < IUABS) exit
          ITUP = KTGEU(ITABS,IUABS)-NTGEUES(ISYM)
          ITUM = KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
          call EXCH(ISYMI,ISYMT,ISYMJ,ISYMU,ITTOT,IUTOT,ERI,SCR)
          if (ITABS /= IUABS) then
            do II=1,NISH(ISYMI)
              IIABS = II+NIES(ISYMI)
              do IJ=1,NISH(ISYMJ)
                IJABS = IJ+NIES(ISYMJ)
                IBUF = II+NORB(ISYMI)*(IJ-1)
                Val = Half*ERI(IBUF)
                if (IIABS >= IJABS) then
                  IIJP = KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                  JWP = ITUP+NASP*(IIJP-1)
                  if (IIABS > IJABS) then
                    GA_Arrays(LWP)%A(JWP) = GA_Arrays(LWP)%A(JWP)+Val
                    IIJM = KIGTJ(IIABS,IJABS)-NIGTJES(ISYM)
                    IWM = ITUM+NASM*(IIJM-1)
                    GA_Arrays(LWM)%A(IWM) = GA_Arrays(LWM)%A(IWM)+Val
                  else
                    GA_Arrays(LWP)%A(JWP) = GA_Arrays(LWP)%A(JWP)+SQ2*Val
                  end if
                else
                  IIJP = KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                  JWP = ITUP+NASP*(IIJP-1)
                  GA_Arrays(LWP)%A(JWP) = GA_Arrays(LWP)%A(JWP)+Val
                  IIJM = KIGTJ(IJABS,IIABS)-NIGTJES(ISYM)
                  IWM = ITUM+NASM*(IIJM-1)
                  GA_Arrays(LWM)%A(IWM) = GA_Arrays(LWM)%A(IWM)-Val
                end if
              end do
            end do
          else
            do II=1,NISH(ISYMI)
              IIABS = II+NIES(ISYMI)
              do IJ=1,NISH(ISYMJ)
                IJABS = IJ+NIES(ISYMJ)
                IBUF = II+NORB(ISYMI)*(IJ-1)
                Val = Quart*ERI(IBUF)
                if (IIABS >= IJABS) then
                  IIJP = KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                  JWP = ITUP+NASP*(IIJP-1)
                  if (IIABS > IJABS) then
                    GA_Arrays(LWP)%A(JWP) = GA_Arrays(LWP)%A(JWP)+Val
                  else
                    GA_Arrays(LWP)%A(JWP) = GA_Arrays(LWP)%A(JWP)+SQ2*Val
                  end if
                else
                  IIJP = KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                  JWP = ITUP+NASP*(IIJP-1)
                  GA_Arrays(LWP)%A(JWP) = GA_Arrays(LWP)%A(JWP)+Val
                end if
              end do
            end do
          end if
        end do
      end do
    end do
  end do
  !  Put WP on disk
  ICASE = 2
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
  ! Put WM on disk
  if (NINM > 0) then
    ICASE = 3
    call MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
  end if
  call Deallocate_GA_Array(LWM)
  call Deallocate_GA_Array(LWP)
end do

end subroutine MKRHSB
