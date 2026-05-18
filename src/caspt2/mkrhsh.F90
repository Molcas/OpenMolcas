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

subroutine MKRHSH(IVEC,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for cases 12 and 13 (BJAI).

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp
use constants, only: half, One, two, three
use SUPERINDEX, only: KAGEB, KIGEJ, KAGTB, KIGTJ
use fake_GA, only: GA_Arrays, Allocate_GA_Array, Deallocate_GA_Array
use caspt2_module, only: NSYM, NAGEB, NIGEJ, NAGTB, NIGTJ, NISH, NIES, NSES, NSSH, NORB, NASH, NAGEBES, NIGEJES, NAGTBES, NIGTJES

implicit none
integer(kind=iwp), intent(in) :: IVEC, nERI1, nERI2, nSCR
real(kind=wp), intent(inout) :: ERI1(nERI1), ERI2(nERI2), SCR(nSCR)
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQI2 = One/SQ2, SQ3 = sqrt(Three)
integer(kind=iwp) ISYM, NASP, NISP, NVP, NASM, NISM, NVM, LVP, ISYMA, ISYMB, ISYMI, ISYMJ, II, IIABS, IJ, IJABS, IA, IAABS, IATOT, &
                  IB, IBABS, IBTOT, IBUF, IVAP, IVIP, IVP, IVAM, IVIM, IVM, LVM, ICASE
real(kind=wp) A, B

do ISYM=1,NSYM
  NASP = NAGEB(ISYM)
  NISP = NIGEJ(ISYM)
  NVP = NASP*NISP
  if (NVP == 0) cycle
  NASM = NAGTB(ISYM)
  NISM = NIGTJ(ISYM)
  NVM = NASM*NISM
  LVP = Allocate_GA_Array(NVP,'WHP')
  if (NVM > 0) LVM = Allocate_GA_Array(NVM,'WHM')
  !  VP(ij,ab)=2*((aibj)+(ajbi))
  ! With new norm., divide by /SQRT(4*(1+Kron(ij))*(1+Kron(ab))
  !   VM(ij,ab)=6*((aibj)-(ajbi))
  ! With new norm., divide by /SQRT(12)
  do ISYMA=1,NSYM
    ISYMB = Mul(ISYMA,ISYM)
    if (ISYMA < ISYMB) cycle
    do ISYMI=1,NSYM
      ISYMJ = Mul(ISYMI,ISYM)
      if (ISYMI < ISYMJ) cycle
      do II=1,NISH(ISYMI)
        IIABS = II+NIES(ISYMI)
        do IJ=1,NISH(ISYMJ)
          IJABS = IJ+NIES(ISYMJ)
          if (IIABS < IJABS) exit
          call EXCH(ISYMA,ISYMI,ISYMB,ISYMJ,II,IJ,ERI1,SCR)
          call EXCH(ISYMA,ISYMJ,ISYMB,ISYMI,IJ,II,ERI2,SCR)
          do IA=1,NSSH(ISYMA)
            IAABS = IA+NSES(ISYMA)
            IATOT = IA+NISH(ISYMA)+NASH(ISYMA)
            do IB=1,NSSH(ISYMB)
              IBABS = IB+NSES(ISYMB)
              if (IAABS < IBABS) exit
              IBTOT = IB+NISH(ISYMB)+NASH(ISYMB)
              IBUF = IATOT+NORB(ISYMA)*(IBTOT-1)
              IVAP = KAGEB(IAABS,IBABS)-NAGEBES(ISYM)
              IVIP = KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
              IVP = IVAP+NAGEB(ISYM)*(IVIP-1)
              A = ERI1(IBUF)+ERI2(IBUF)
              if (IIABS /= IJABS) then
                if (IAABS /= IBABS) then
                  GA_Arrays(LVP)%A(IVP) = A
                  IVAM = KAGTB(IAABS,IBABS)-NAGTBES(ISYM)
                  IVIM = KIGTJ(IIABS,IJABS)-NIGTJES(ISYM)
                  IVM = IVAM+NAGTB(ISYM)*(IVIM-1)
                  B = ERI1(IBUF)-ERI2(IBUF)
                  GA_Arrays(LVM)%A(IVM) = SQ3*B
                else
                  GA_Arrays(LVP)%A(IVP) = SQI2*A
                end if
              else
                if (IAABS /= IBABS) then
                  GA_Arrays(LVP)%A(IVP) = SQI2*A
                else
                  GA_Arrays(LVP)%A(IVP) = Half*A
                end if
              end if
            end do
          end do
        end do
      end do
    end do
  end do

  ICASE = 12
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LVP)
  call Deallocate_GA_Array(LVP)
  if (NVM > 0) then
    ICASE = 13
    call MKRHS_SAVE(ICASE,ISYM,IVEC,LVM)
    call Deallocate_GA_Array(LVM)
  end if
end do

end subroutine MKRHSH
