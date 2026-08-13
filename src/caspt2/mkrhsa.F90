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

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use SUPERINDEX, only: KTUV
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array, GA_Arrays
use general_data, only: NACTEL, NASH
use caspt2_module, only: NAES, NINDEP, NISH, NISH, NORB, NSYM, NTUV, NTUVES
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, NFIMO, nERI, nSCR
real(kind=wp), intent(inout) :: FIMO(NFIMO), ERI(nERI), SCR(nSCR)
integer(kind=iwp) :: IBUF, ICASE, II, ISYM, ISYMT, ISYMU, ISYMUV, ISYMV, IT, ITABS, ITTOT, IU, IUABS, IUTOT, IV, IVABS, IVTOT, IW, &
                     IW1, IW2, LW, NAS, NFIMOES, NFNXT, NI, NIS, NV
real(kind=wp) :: FTI, ONEADD, WTUVI

NFNXT = 0
do ISYM=1,NSYM
  NFIMOES = NFNXT
  NFNXT = NFNXT+nTri_Elem(NORB(ISYM))
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
            FTI = FIMO(NFIMOES+iTri(ITTOT,II))
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
