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

subroutine MKRHSE(IVEC,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for cases 6 and 7 (VJAI).

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp
use constants, only: half, One, two, three
use SUPERINDEX, only: KIGEJ, KIGTJ
use fake_GA, only: GA_Arrays, Allocate_GA_Array, Deallocate_GA_Array
use caspt2_module, only: NSYM, NINDEP, NISUP, NASH, NISH, NSSH, NORB, NIGEJ, NIES, NIGEJES, NIGTJES, NIGTJ

implicit none
integer(kind=iwp), intent(in) :: IVEC, nERI1, nERI2, nSCR
real(kind=wp), intent(inout) :: ERI1(nERI1), ERI2(nERI2), SCR(nSCR)
integer(kind=iwp) IOFF1(8), IOFF2(8)
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQI2 = One/SQ2, SQ3 = sqrt(Three), SQ32 = SQ3*SQI2
integer(kind=iwp) ISYM, IO1, IO2, NAS, NISP, NISM, NVP, NVM, LWP, LWM, ISYMA, ISYMI, IT, ITTOT, II, IA, IGEJ, IGTJ, IIABS, IATOT, &
                  IBUF, IWA, IWIP, JWP, IJ, IJABS, ISYMIJ, ISYMJ, IWIM, IWM, ICASE
real(kind=wp) A, B

do ISYM=1,NSYM
  if (NINDEP(ISYM,6)+NINDEP(ISYM,7) == 0) cycle
  ! Set up offset table:
  IO1 = 0
  IO2 = 0
  do ISYMA=1,NSYM
    IOFF1(ISYMA) = IO1
    IOFF2(ISYMA) = IO2
    ISYMIJ = Mul(ISYMA,ISYM)
    IO1 = IO1+NSSH(ISYMA)*NIGEJ(ISYMIJ)
    IO2 = IO2+NSSH(ISYMA)*NIGTJ(ISYMIJ)
  end do
  ! Allocate W with parts WP,WM
  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,6)
  NISM = NISUP(ISYM,7)
  NVP = NAS*NISP
  if (NVP == 0) cycle
  NVM = NAS*NISM
  LWP = Allocate_GA_Array(NVP,'WEP')
  LWM = Allocate_GA_Array(NVM,'WEM')
  ! Let W(t,i,j,a)=(aitj)
  !   WP(t,ij,a)=  (W(t,i,j,a)+W(t,j,i,a))
  ! With new normalisation, divide by /SQRT(2+2*Kron(ij))
  !   WM(t,ij,a)=3*(W(t,i,j,a)-W(t,j,i,a))
  ! With new normalisation, divide by /SQRT(6)
  do ISYMA=1,NSYM
    ISYMIJ = Mul(ISYMA,ISYM)
    do ISYMI=1,NSYM
      ISYMJ = Mul(ISYMI,ISYMIJ)
      if (ISYMI < ISYMJ) cycle
      do II=1,NISH(ISYMI)
        IIABS = II+NIES(ISYMI)
        do IJ=1,NISH(ISYMJ)
          IJABS = IJ+NIES(ISYMJ)
          if (IIABS < IJABS) exit
          call EXCH(ISYMA,ISYMI,ISYM,ISYMJ,II,IJ,ERI1,SCR)
          call EXCH(ISYMA,ISYMJ,ISYM,ISYMI,IJ,II,ERI2,SCR)
          IGEJ = KIGEJ(IIABS,IJABS)-NIGEJES(ISYMIJ)
          IGTJ = KIGTJ(IIABS,IJABS)-NIGTJES(ISYMIJ)
          do IA=1,NSSH(ISYMA)
            IATOT = IA+NISH(ISYMA)+NASH(ISYMA)
            do IT=1,NASH(ISYM)
              ITTOT = IT+NISH(ISYM)
              IBUF = IATOT+NORB(ISYMA)*(ITTOT-1)
              A = ERI1(IBUF)+ERI2(IBUF)
              IWA = IT
              IWIP = IA+NSSH(ISYMA)*(IGEJ-1)+IOFF1(ISYMA)
              JWP = IWA+NAS*(IWIP-1)
              if (IIABS > IJABS) then
                GA_Arrays(LWP)%A(JWP) = SQI2*A
                B = ERI1(IBUF)-ERI2(IBUF)
                IWIM = IA+NSSH(ISYMA)*(IGTJ-1)+IOFF2(ISYMA)
                IWM = IWA+NAS*(IWIM-1)
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
  ICASE = 6
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
  if (NVM > 0) then
    ICASE = 7
    call MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
  end if
  call Deallocate_GA_Array(LWP)
  call Deallocate_GA_Array(LWM)
end do

end subroutine MKRHSE
