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

subroutine MKRHSC(IVEC,FIMO,NFIMO,ERI,nERI,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV for case 4 (ATVX).

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp
use SUPERINDEX, only: KTUV
use fake_GA, only: GA_Arrays, Allocate_GA_Array, Deallocate_GA_Array
use caspt2_module, only: NSYM, NORB, NINDEP, NTUV, NSSH, NASH, NISH, NAES, NSSH, NTUVES, NASHT, NACTEL

implicit none
integer(kind=iwp), intent(in) :: IVEC, NFIMO, nERI, nSCR
real(kind=wp), intent(inout) :: FIMO(NFIMO), ERI(nERI), SCR(nSCR)
integer(kind=iwp) NFNXT, ISYM, NFIMOES, NAS, NIS, NV, LW, ISYMT, ISYMUV, ISYMU, ISYMV, IU, IUTOT, IUABS, IV, IVTOT, IVABS, IA, &
                  IATOT, IT, ITTOT, ITABS, IW1, IW2, IW, IBUF, IFIMO, IYABS, IYYW, IYYWA, ICASE
real(kind=wp) SUM, ONEADD

NFNXT = 0
do ISYM=1,NSYM
  NFIMOES = NFNXT
  NFNXT = NFNXT+(NORB(ISYM)*(NORB(ISYM)+1))/2
  if (NINDEP(ISYM,4) == 0) cycle
  NAS = NTUV(ISYM)
  NIS = NSSH(ISYM)
  NV = NAS*NIS
  if (NV == 0) cycle

  !   Allocate W. Put in W(tuv,a)=(at,uv) +
  !             (FIMO(a,t)-sum(y)(ay,yt))*delta(u,v)/NACTEL.
  ! First, just the two-electron integrals. Later, add correction.

  LW = Allocate_GA_Array(NV,'WC')
  do ISYMT=1,NSYM
    ISYMUV = Mul(ISYMT,ISYM)
    do ISYMU=1,NSYM
      ISYMV = Mul(ISYMU,ISYMUV)
      do IU=1,NASH(ISYMU)
        IUTOT = IU+NISH(ISYMU)
        IUABS = IU+NAES(ISYMU)
        do IV=1,NASH(ISYMV)
          IVTOT = IV+NISH(ISYMV)
          IVABS = IV+NAES(ISYMV)
          call COUL(ISYM,ISYMT,ISYMU,ISYMV,IUTOT,IVTOT,ERI,SCR)
          do IA=1,NSSH(ISYM)
            IATOT = IA+NISH(ISYM)+NASH(ISYM)
            do IT=1,NASH(ISYMT)
              ITTOT = IT+NISH(ISYMT)
              ITABS = IT+NAES(ISYMT)
              IW1 = KTUV(ITABS,IUABS,IVABS)-NTUVES(ISYM)
              IW2 = IA
              IW = IW1+NAS*(IW2-1)
              IBUF = IATOT+NORB(ISYM)*(ITTOT-1)
              GA_Arrays(LW)%A(IW) = ERI(IBUF)
            end do
          end do
        end do
      end do
    end do
  end do

  do IT=1,NASH(ISYM)
    ITTOT = IT+NISH(ISYM)
    ITABS = IT+NAES(ISYM)
    do IA=1,NSSH(ISYM)
      IATOT = IA+NISH(ISYM)+NASH(ISYM)
      IFIMO = NFIMOES+(IATOT*(IATOT-1))/2+ITTOT
      SUM = FIMO(IFIMO)
      do IYABS=1,NASHT
        IYYW = KTUV(IYABS,IYABS,ITABS)-NTUVES(ISYM)
        IYYWA = IYYW+NAS*(IA-1)
        SUM = SUM-GA_Arrays(LW)%A(IYYWA)
      end do
      ONEADD = SUM/dble(max(1,NACTEL))
      do ISYMU=1,NSYM
        do IU=1,NASH(ISYMU)
          IUABS = IU+NAES(ISYMU)
          IW1 = KTUV(ITABS,IUABS,IUABS)-NTUVES(ISYM)
          IW2 = IA
          IW = IW1+NAS*(IW2-1)
          GA_Arrays(LW)%A(IW) = GA_Arrays(LW)%A(IW)+ONEADD
        end do
      end do
    end do
  end do

  ! Put W on disk
  ICASE = 4
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LW)

  call Deallocate_GA_Array(LW)
end do

end subroutine MKRHSC
