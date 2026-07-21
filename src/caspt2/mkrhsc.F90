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

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use SUPERINDEX, only: KTUV
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array, GA_Arrays
use general_data, only: NACTEL, NASH
use caspt2_module, only: NAES, NASHT, NINDEP, NISH, NORB, NSSH, NSSH, NSYM, NTUV, NTUVES
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, NFIMO, nERI, nSCR
real(kind=wp), intent(inout) :: FIMO(NFIMO), ERI(nERI), SCR(nSCR)

integer(kind=iwp) :: IA, IATOT, IBUF, ICASE, IFIMO, ISYM, ISYMT, ISYMU, ISYMUV, ISYMV, IT, ITABS, ITTOT, IU, IUABS, IUTOT, IV, &
                     IVABS, IVTOT, IW, IW1, IW2, IYABS, IYYW, IYYWA, LW, NAS, NFIMOES, NFNXT, NIS, NV
real(kind=wp) :: ONEADD, rSUM

NFNXT = 0
do ISYM=1,NSYM
  NFIMOES = NFNXT
  NFNXT = NFNXT+nTri_Elem(NORB(ISYM))
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
      IFIMO = NFIMOES+iTri(IATOT,ITTOT)
      rSUM = FIMO(IFIMO)
      do IYABS=1,NASHT
        IYYW = KTUV(IYABS,IYABS,ITABS)-NTUVES(ISYM)
        IYYWA = IYYW+NAS*(IA-1)
        rSUM = rSUM-GA_Arrays(LW)%A(IYYWA)
      end do
      ONEADD = rSUM/real(max(1,NACTEL),kind=wp)
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
