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

subroutine MKRHSA(IVEC,FIMO,NFIMO,ERI,nERI,SCR,nSCR)
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for case 1 (VJTU).

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp
use constants, only: Zero
use SUPERINDEX, only: KTUV
use fake_GA, only: GA_Arrays, Allocate_GA_Array, Deallocate_GA_Array
use caspt2_module, only: NSYM, NORB, NINDEP, NTUV, NISH, NASH, NISH, NAES, NTUVES, NACTEL

implicit none
integer(kind=iwp), intent(in) :: IVEC, NFIMO, nERI, nSCR
real(kind=wp), intent(inout) :: FIMO(NFIMO), ERI(nERI), SCR(nSCR)
integer(kind=iwp) NFNXT, ISYM, NFIMOES, NAS, NIS, NV, NI, LW, ISYMT, ISYMUV, ISYMU, ISYMV, IT, ITTOT, ITABS, II, IU, IUTOT, IUABS, &
                  IV, IVTOT, IVABS, IW1, IW2, IW, IBUF, ICASE
real(kind=wp) FTI, ONEADD, WTUVI

NFNXT = 0
do ISYM=1,NSYM
  NFIMOES = NFNXT
  NFNXT = NFNXT+(NORB(ISYM)*(NORB(ISYM)+1))/2
  if (NINDEP(ISYM,1) == 0) cycle
  NAS = NTUV(ISYM)
  NIS = NISH(ISYM)
  NV = NAS*NIS
  if (NV == 0) cycle
  ! Set up a matrix FWI(w,i)=FIMO(wi)
  NI = NISH(ISYM)

  ! Compute W(tuv,i)=(ti,uv) + FIMO(t,i)*delta(u,v)/NACTEL
  LW = Allocate_GA_Array(NV,'WA')
  do ISYMT=1,NSYM
    ISYMUV = Mul(ISYMT,ISYM)
    do ISYMU=1,NSYM
      ISYMV = Mul(ISYMU,ISYMUV)
      do IT=1,NASH(ISYMT)
        ITTOT = IT+NISH(ISYMT)
        ITABS = IT+NAES(ISYMT)
        do II=1,NI
          call COUL(ISYMU,ISYMV,ISYMT,ISYM,ITTOT,II,ERI,SCR)
          ONEADD = Zero
          if (ISYMT == ISYM) then
            FTI = FIMO(NFIMOES+(ITTOT*(ITTOT-1))/2+II)
            ONEADD = FTI/real(max(1,NACTEL),kind=wp)
          end if
          do IU=1,NASH(ISYMU)
            IUTOT = IU+NISH(ISYMU)
            IUABS = IU+NAES(ISYMU)
            do IV=1,NASH(ISYMV)
              IVTOT = IV+NISH(ISYMV)
              IVABS = IV+NAES(ISYMV)
              IW1 = KTUV(ITABS,IUABS,IVABS)-NTUVES(ISYM)
              IW2 = II
              IW = IW1+NAS*(IW2-1)
              IBUF = IUTOT+NORB(ISYMU)*(IVTOT-1)
              WTUVI = ERI(IBUF)
              if (IVABS == IUABS) WTUVI = WTUVI+ONEADD
              GA_Arrays(LW)%A(IW) = WTUVI
            end do
          end do
        end do
      end do
    end do
  end do
  ! Put W on disk:
  ICASE = 1
  call MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
  call Deallocate_GA_Array(LW)
end do

end subroutine MKRHSA
