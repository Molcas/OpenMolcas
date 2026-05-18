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

subroutine MKRHSG(IVEC,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for cases 10 and 11 (BJAT).

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp
use constants, only: half, One, two, three
use SUPERINDEX, only: KAGEB, KAGTB
use fake_GA, only: GA_Arrays, Allocate_GA_Array, Deallocate_GA_Array
use caspt2_module, only: NSYM, NINDEP, NISH, NAGEB, NAGTB, NASH, NISUP, NISUP, NSSH, NSES, NORB, NAGEBES, NAGTBES

implicit none
integer(kind=iwp), intent(in) :: IVEC, nERI1, nERI2, nSCR
real(kind=wp), intent(inout) :: ERI1(nERI1), ERI2(nERI2), SCR(nSCR)
integer(kind=iwp) IOFF1(8), IOFF2(8)
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQI2 = One/SQ2, SQ3 = sqrt(Three), SQ32 = SQ3*SQI2
integer(kind=iwp) ISYM, IO1, IO2, NAS, NISP, NISM, NVP, NVM, LWP, LWM, ISYMA, ISYMB, ISYMAB, ISYMI, IT, ITTOT, II, IA, IAABS, &
                  IATOT, IB, IBABS, IBTOT, IBUF, IWA, IAGEB, IWIP, JWP, IAGTB, IWIM, IWM, ICASE
real(kind=wp) A, B

do ISYM=1,NSYM
  if (NINDEP(ISYM,10)+NINDEP(ISYM,11) == 0) cycle
  ! Set up offset table:
  IO1 = 0
  IO2 = 0
  do ISYMI=1,NSYM
    IOFF1(ISYMI) = IO1
    IOFF2(ISYMI) = IO2
    ISYMAB = Mul(ISYMI,ISYM)
    IO1 = IO1+NISH(ISYMI)*NAGEB(ISYMAB)
    IO2 = IO2+NISH(ISYMI)*NAGTB(ISYMAB)
  end do
  ! Allocate W with parts WP,WM
  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,10)
  NISM = NISUP(ISYM,11)
  NVP = NAS*NISP
  if (NVP == 0) cycle
  NVM = NAS*NISM
  LWP = Allocate_GA_Array(NVP,'WGP')
  LWM = Allocate_GA_Array(NVM,'WGM')
  ! Let  W(t,i,a,b)=(atbi)
  !   WP(t,i,ab)=  (W(t,i,a,b)+W(t,i,b,a))
  ! With new normalisation, divide by /SQRT(2+2*Kron(ab))
  !   WM(t,i,ab)=3*(W(t,i,a,b)-W(t,i,b,a))
  ! With new normalisation, divide by /SQRT(6)
  do ISYMA=1,NSYM
    do ISYMB=1,ISYMA
      ISYMAB = Mul(ISYMA,ISYMB)
      ISYMI = Mul(ISYMAB,ISYM)
      do IT=1,NASH(ISYM)
        ITTOT = IT+NISH(ISYM)
        do II=1,NISH(ISYMI)
          call EXCH(ISYMA,ISYM,ISYMB,ISYMI,ITTOT,II,ERI1,SCR)
          call EXCH(ISYMA,ISYMI,ISYMB,ISYM,II,ITTOT,ERI2,SCR)
          do IA=1,NSSH(ISYMA)
            IAABS = IA+NSES(ISYMA)
            IATOT = IA+NISH(ISYMA)+NASH(ISYMA)
            do IB=1,NSSH(ISYMB)
              IBABS = IB+NSES(ISYMB)
              IBTOT = IB+NISH(ISYMB)+NASH(ISYMB)
              if (IAABS < IBABS) exit
              IBUF = IATOT+NORB(ISYMA)*(IBTOT-1)
              IWA = IT
              IAGEB = KAGEB(IAABS,IBABS)-NAGEBES(ISYMAB)
              IWIP = II+NISH(ISYMI)*(IAGEB-1)+IOFF1(ISYMI)
              JWP = IWA+NAS*(IWIP-1)
              A = ERI1(IBUF)+ERI2(IBUF)
              if (IAABS /= IBABS) then
                GA_Arrays(LWP)%A(JWP) = SQI2*A
                IAGTB = KAGTB(IAABS,IBABS)-NAGTBES(ISYMAB)
                IWIM = II+NISH(ISYMI)*(IAGTB-1)+IOFF2(ISYMI)
                IWM = IWA+NAS*(IWIM-1)
                B = ERI1(IBUF)-ERI2(IBUF)
                GA_Arrays(LWM)%A(IWM) = SQ32*B
              else
                GA_Arrays(LWP)%A(JWP) = half*A
              end if
            end do
          end do
        end do
      end do
    end do
  end do
  ! Put WP and WM on disk.
  ICASE = 10
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
  if (NVM > 0) then
    ICASE = 11
    call MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
  end if
  call Deallocate_GA_Array(LWP)
  call Deallocate_GA_Array(LWM)
end do

end subroutine MKRHSG
