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
! Copyright (C) Steven Vancoillie                                      *
!***********************************************************************

subroutine RHSOD_G_NOSYM(IVEC)

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp, u6
use constants, only: Half, OneHalf
use SUPERINDEX, only: MAGEB, MAREL, MAGTB
use CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG
use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
use fake_GA, only: GA_Arrays
#endif
use caspt2_module, only: NSYM, NASUP, NISUP, NISH, NAGEB, NAGEBES, NSSH, NAGTB, NAGTBES

implicit none
integer(iwp), intent(in) :: IVEC
integer(iwp) IOBRA(8,8), IOKET(8,8)
real(kind=wp), allocatable :: BRABUF(:), KETBUF(:)
real(kind=wp), parameter :: SQRTH = sqrt(Half), SQRTA = sqrt(OneHalf)
real(kind=wp) AVCJ, CVAJ, GM, GP, SCL
integer(iwp) IA, IAABS, NAS, NIS, lg_W, IASTA, IAEND, IISTA, IIEND, MW, IAGEC, IAGECTOT, IAGTC, IAGTCTOT, IAJ, IAV, IC, ICABS, &
             iCASE, ICJ, ICV, IDX, IJ, IJAGEC, IJAGECEND, IJAGECSTA, IJAGTC, IJAGTCEND, IJAGTCSTA, IOFF, IOFFAJ, IOFFAV, IOFFCJ, &
             IOFFCV, ISYA, ISYAC, ISYC, ISYJ, ISYM, ISYV, IV, NAC, NBRABUF, NJ, NKETBUF, NV, NW
real(kind=wp), external :: DDOT_
!logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

if (iPrGlb >= DEBUG) write(u6,*) 'RHS on demand: case G'

!***********************************************************************
! Case G (10,11):
! GP(v,jac)=((av,cj)+(cv,aj))/SQRT(2+2*Kron(a,b))
! GM(v,jac)=((av,cj)-(cv,aj))*SQRT(3/2)
!***********************************************************************

! -SVC- Case G is slightly special, in that the inactive superindices are
! so large, that it is suboptimal to have a direct translation table for
! them. Instead, the code loops over symmetry blocks of J-AC and figures
! out if the indices on the processor fall within a block or not. Within
! a J-AC symmetry block, NJ(ISYJ) and NAGEB(ISYAC) are known, so they can
! be determined by integer division. This could be optimized by combining
! it with loop peeling (on the todo list?).

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
call CHOVEC_SIZE(3,NBRABUF,IOBRA)
call CHOVEC_SIZE(4,NKETBUF,IOKET)

call mma_allocate(BRABUF,NBRABUF,LABEL='BRABUF')
call mma_allocate(KETBUF,NKETBUF,LABEL='KETBUF')

call CHOVEC_READ(3,BRABUF,NBRABUF)
call CHOVEC_READ(4,KETBUF,NKETBUF)

iCASE = 10
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
do ISYM=1,NSYM

  NAS = NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE)
  NW = NAS*NIS

  if (NW == 0) cycle

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
  NW = NAS*(IIEND-IISTA+1)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  ! find start and end block
  IOFF = 0
  do ISYJ=1,NSYM
    ISYAC = Mul(ISYJ,ISYM)
    ISYV = ISYM

    NJ = NISH(ISYJ)
    NAC = NAGEB(ISYAC)
    ! what is start/end in this block?
    IJAGECSTA = max(IISTA-IOFF,1)
    IJAGECEND = min(IIEND-IOFF,NJ*NAC)

    do IJAGEC=IJAGECSTA,IJAGECEND
      IAGEC = (IJAGEC-1)/NJ+1
      IJ = IJAGEC-NJ*(IAGEC-1)
      IAGECTOT = IAGEC+NAGEBES(ISYAC)
      IAABS = MAGEB(1,IAGECTOT)
      ICABS = MAGEB(2,IAGECTOT)
      IA = MAREL(1,IAABS)
      ISYA = MAREL(2,IAABS)
      IC = MAREL(1,ICABS)
      ISYC = MAREL(2,ICABS)
      do IV=IASTA,IAEND ! these are always all elements
        ! compute integrals (ajvl) and (alvj)
        NV = NVTOT_CHOSYM(Mul(ISYA,ISYV))
        IAV = IA-1+NSSH(ISYA)*(IV-1)
        ICJ = IC-1+NSSH(ISYC)*(IJ-1)
        IOFFAV = 1+IOBRA(ISYA,ISYV)+NV*IAV
        IOFFCJ = 1+IOKET(ISYC,ISYJ)+NV*ICJ
        AVCJ = DDOT_(NV,BRABUF(IOFFAV),1,KETBUF(IOFFCJ),1)

        NV = NVTOT_CHOSYM(Mul(ISYC,ISYV))
        ICV = IC-1+NSSH(ISYC)*(IV-1)
        IAJ = IA-1+NSSH(ISYA)*(IJ-1)
        IOFFCV = 1+IOBRA(ISYC,ISYV)+NV*ICV
        IOFFAJ = 1+IOKET(ISYA,ISYJ)+NV*IAJ
        CVAJ = DDOT_(NV,BRABUF(IOFFCV),1,KETBUF(IOFFAJ),1)

        ! GP(v,jac)=((av,cj)+(cv,aj))/SQRT(2+2*Kron(a,b))
        if (IAABS == ICABS) then
          SCL = Half
        else
          SCL = SQRTH
        end if
        GP = SCL*(AVCJ+CVAJ)
        ! write element EP
        IDX = IV+NAS*(IJAGEC+IOFF-IISTA)
#       ifdef _MOLCAS_MPP_
        DBL_MB(MW+IDX-1) = GP
#       else
        GA_Arrays(lg_w)%A(IDX) = GP
#       endif
      end do
    end do

    IOFF = IOFF+NJ*NAC
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

iCASE = 11
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
do ISYM=1,NSYM

  NAS = NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE)
  NW = NAS*NIS

  if (NW == 0) cycle

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
  NW = NAS*(IIEND-IISTA+1)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  ! find start and end block
  IOFF = 0
  do ISYJ=1,NSYM
    ISYAC = Mul(ISYJ,ISYM)
    ISYV = ISYM

    NJ = NISH(ISYJ)
    NAC = NAGTB(ISYAC)
    ! what is start/end in this block?
    IJAGTCSTA = max(IISTA-IOFF,1)
    IJAGTCEND = min(IIEND-IOFF,NJ*NAC)

    do IJAGTC=IJAGTCSTA,IJAGTCEND
      IAGTC = (IJAGTC-1)/NJ+1
      IJ = IJAGTC-NJ*(IAGTC-1)
      IAGTCTOT = IAGTC+NAGTBES(ISYAC)
      IAABS = MAGTB(1,IAGTCTOT)
      ICABS = MAGTB(2,IAGTCTOT)
      IA = MAREL(1,IAABS)
      ISYA = MAREL(2,IAABS)
      IC = MAREL(1,ICABS)
      ISYC = MAREL(2,ICABS)
      do IV=IASTA,IAEND ! these are always all elements
        ! compute integrals (ajvl) and (alvj)
        NV = NVTOT_CHOSYM(Mul(ISYA,ISYV))
        IAV = IA-1+NSSH(ISYA)*(IV-1)
        ICJ = IC-1+NSSH(ISYC)*(IJ-1)
        IOFFAV = 1+IOBRA(ISYA,ISYV)+NV*IAV
        IOFFCJ = 1+IOKET(ISYC,ISYJ)+NV*ICJ
        AVCJ = DDOT_(NV,BRABUF(IOFFAV),1,KETBUF(IOFFCJ),1)

        NV = NVTOT_CHOSYM(Mul(ISYC,ISYV))
        ICV = IC-1+NSSH(ISYC)*(IV-1)
        IAJ = IA-1+NSSH(ISYA)*(IJ-1)
        IOFFCV = 1+IOBRA(ISYC,ISYV)+NV*ICV
        IOFFAJ = 1+IOKET(ISYA,ISYJ)+NV*IAJ
        CVAJ = DDOT_(NV,BRABUF(IOFFCV),1,KETBUF(IOFFAJ),1)

        ! GM(v,jac)=((av,cj)-(cv,aj))*SQRT(3/2)
        GM = SQRTA*(AVCJ-CVAJ)
        ! write element GM
        IDX = IV+NAS*(IJAGTC+IOFF-IISTA)
#       ifdef _MOLCAS_MPP_
        DBL_MB(MW+IDX-1) = GM
#       else
        GA_Arrays(lg_w)%A(IDX) = GM
#       endif
      end do
    end do

    IOFF = IOFF+NJ*NAC
  end do
  !*********************************************************************

  call RHS_Release_Update(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

call mma_deallocate(BRABUF)
call mma_deallocate(KETBUF)

end subroutine RHSOD_G_NOSYM
