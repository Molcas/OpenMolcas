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

subroutine MKRHSD(IVEC,FIMO,NFIMO,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for case 5, AIVX.

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp
use constants, only: Zero
use SUPERINDEX, only: KTU
use fake_GA, only: GA_Arrays, Allocate_GA_Array, Deallocate_GA_Array
use caspt2_module, only: NSYM, NINDEP, NSSH, NISH, NTU, NISUP, NACTEL, NORB, NAES, NASH, NTUES

implicit none
integer(kind=iwp), intent(in) :: IVEC, NFIMO, nERI1, nERI2, nSCR
real(kind=wp), intent(inout) :: FIMO(NFIMO)
real(kind=wp), intent(inout) :: ERI1(nERI1), ERI2(nERI2), SCR(nSCR)
integer(kind=iwp) IOFF(8)
integer(kind=iwp) ISYM, IO, ISYMI, NAS1, NAS, NIS, NV, LW, NFSUM, NFIMOES, ISYMA, ISYMU, ISYMT, II, IU, IUABS, IUTOT, IA, IATOT, &
                  IT, ITABS, ITTOT, IWA, IWI, IW1, IW2, IBUF1, IBUF2, ICASE
real(kind=wp) ONEADD, FAI, WAITU

do ISYM=1,NSYM
  if (NINDEP(ISYM,5) == 0) cycle
  ! Set up offset table:
  IO = 0
  do ISYMA=1,NSYM
    IOFF(ISYMA) = IO
    ISYMI = Mul(ISYMA,ISYM)
    IO = IO+NSSH(ISYMA)*NISH(ISYMI)
  end do
  ! Allocate W; W subdivided into W1,W2.
  NAS1 = NTU(ISYM)
  NAS = 2*NAS1
  NIS = NISUP(ISYM,5)
  NV = NAS*NIS
  if (NV == 0) cycle
  ! Compute W1(tu,ai)=(ai,tu) + FIMO(a,i)*delta(t,u)/NACTEL
  ! Compute W2(tu,ai)=(ti,au)
  LW = Allocate_GA_Array(NV,'WD')
  NFSUM = 0
  do ISYMI=1,NSYM
    NFIMOES = NFSUM
    NFSUM = NFSUM+(NORB(ISYMI)*(NORB(ISYMI)+1))/2
    ISYMA = Mul(ISYMI,ISYM)
    do ISYMU=1,NSYM
      ISYMT = Mul(ISYMU,ISYM)
      do II=1,NISH(ISYMI)
        do IU=1,NASH(ISYMU)
          IUABS = IU+NAES(ISYMU)
          IUTOT = IU+NISH(ISYMU)
          call EXCH(ISYMA,ISYMI,ISYMT,ISYMU,II,IUTOT,ERI1,SCR)
          call EXCH(ISYMT,ISYMI,ISYMA,ISYMU,II,IUTOT,ERI2,SCR)
          do IA=1,NSSH(ISYMA)
            IATOT = IA+NISH(ISYMA)+NASH(ISYMA)
            ONEADD = Zero
            if (ISYM == 1) then
              FAI = FIMO(NFIMOES+(IATOT*(IATOT-1))/2+II)
              ONEADD = FAI/dble(max(1,NACTEL))
            end if
            do IT=1,NASH(ISYMT)
              ITABS = IT+NAES(ISYMT)
              ITTOT = IT+NISH(ISYMT)
              IWA = KTU(ITABS,IUABS)-NTUES(ISYM)
              IWI = II+NISH(ISYMI)*(IA-1)+IOFF(ISYMA)
              IW1 = IWA+NAS*(IWI-1)
              IW2 = IW1+NAS1
              IBUF1 = IATOT+NORB(ISYMA)*(ITTOT-1)
              IBUF2 = ITTOT+NORB(ISYMT)*(IATOT-1)
              WAITU = ERI1(IBUF1)
              if (ITABS == IUABS) WAITU = WAITU+ONEADD
              GA_Arrays(LW)%A(IW1) = WAITU
              GA_Arrays(LW)%A(IW2) = ERI2(IBUF2)
            end do
          end do
        end do
      end do
    end do
  end do
  ! Put W on disk.
  ICASE = 5
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
  call Deallocate_GA_Array(LW)
end do

end subroutine MKRHSD
