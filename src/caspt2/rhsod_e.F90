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

subroutine RHSOD_E(IVEC)

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp, u6
use constants, only: Half, OneHalf
use SUPERINDEX, only: MIGEJ, MIREL, MIGTJ
use CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG
use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
use fake_GA, only: GA_Arrays
#endif
use caspt2_module, only: NSYM, NASUP, NISUP, NSSH, NIGEJ, NIGEJES, NASH, NIGTJ, NIGTJES

implicit none
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) IOBRA(8,8), IOKET(8,8)
real(kind=wp), allocatable :: BRABUF(:), KETBUF(:)
real(kind=wp), parameter :: SQRTH = sqrt(Half), SQRTA = sqrt(OneHalf)
real(kind=wp) AJVL, ALVJ, EM, EP, SCL
integer(kind=iwp) NAS, NIS, lg_W, IASTA, IAEND, IISTA, IIEND, MW, IA, IAJ, IAJGEL, IAJGELEND, IAJGELSTA, IAJGTL, IAJGTLEND, &
                  IAJGTLSTA, IAL, iCASE, IDX, IJ, IJABS, IJGEL, IJGELTOT, IJGTL, IJGTLTOT, IL, ILABS, IOFF, IOFFAJ, IOFFAL, &
                  IOFFVJ, IOFFVL, ISYA, ISYJ, ISYJL, ISYL, ISYM, ISYV, IV, IVJ, IVL, NA, NBRABUF, NJL, NKETBUF, NV, NW
real(kind=wp), external :: DDot_
!logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

if (iPrGlb >= DEBUG) write(u6,*) 'RHS on demand: case E'

!***********************************************************************
! Case E (6,7):
! EP(v,ajl)=((aj,vl)+(al,vj))/SQRT(2+2*Kron(j,l))
! EM(v,ajl)=((aj,vl)-(al,vj))*SQRT(3/2)
!***********************************************************************

! -SVC- Case E is slightly special, in that the inactive superindices are
! so large, that it is suboptimal to have a direct translation table for
! them. Instead, the code loops over symmetry blocks of A-JL and figures
! out if the indices on the processor fall within a block or not. Within
! a A-JL symmetry block, NA(ISYA) and NIGEJ(ISYJL) are known, so they can
! be determined by integer(kind=iwp) division. This could be optimized by combining
! it with loop peeling (on the todo list?).

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
call CHOVEC_SIZE(4,NBRABUF,IOBRA)
call CHOVEC_SIZE(1,NKETBUF,IOKET)

call mma_allocate(BRABUF,NBRABUF,LABEL='BRABUF')
call mma_allocate(KETBUF,NKETBUF,LABEL='KETBUF')

call CHOVEC_READ(4,BRABUF,NBRABUF)
call CHOVEC_READ(1,KETBUF,NKETBUF)

iCASE = 6
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
  do ISYA=1,NSYM
    ISYJL = Mul(ISYA,ISYM)
    ISYV = ISYM

    NA = NSSH(ISYA)
    NJL = NIGEJ(ISYJL)
    ! what is start/end in this block?
    IAJGELSTA = max(IISTA-IOFF,1)
    IAJGELEND = min(IIEND-IOFF,NA*NJL)

    do IAJGEL=IAJGELSTA,IAJGELEND
      IJGEL = (IAJGEL-1)/NA+1
      IA = IAJGEL-NA*(IJGEL-1)
      IJGELTOT = IJGEL+NIGEJES(ISYJL)
      IJABS = MIGEJ(1,IJGELTOT)
      ILABS = MIGEJ(2,IJGELTOT)
      IJ = MIREL(1,IJABS)
      ISYJ = MIREL(2,IJABS)
      IL = MIREL(1,ILABS)
      ISYL = MIREL(2,ILABS)
      do IV=IASTA,IAEND ! these are always all elements
        ! compute integrals (ajvl) and (alvj)
        NV = NVTOT_CHOSYM(Mul(ISYA,ISYJ))
        IAJ = IA-1+NSSH(ISYA)*(IJ-1)
        IVL = IV-1+NASH(ISYV)*(IL-1)
        IOFFAJ = 1+IOBRA(ISYA,ISYJ)+NV*IAJ
        IOFFVL = 1+IOKET(ISYV,ISYL)+NV*IVL
        AJVL = DDOT_(NV,BRABUF(IOFFAJ),1,KETBUF(IOFFVL),1)

        NV = NVTOT_CHOSYM(Mul(ISYA,ISYL))
        IAL = IA-1+NSSH(ISYA)*(IL-1)
        IVJ = IV-1+NASH(ISYV)*(IJ-1)
        IOFFAL = 1+IOBRA(ISYA,ISYL)+NV*IAL
        IOFFVJ = 1+IOKET(ISYV,ISYJ)+NV*IVJ
        ALVJ = DDOT_(NV,BRABUF(IOFFAL),1,KETBUF(IOFFVJ),1)

        ! EP(v,ajl)=((aj,vl)+(al,vj))/SQRT(2+2*Kron(j,l))
        if (ILABS == IJABS) then
          SCL = Half
        else
          SCL = SQRTH
        end if
        EP = SCL*(AJVL+ALVJ)
        ! write element EP
        IDX = IV+NAS*(IAJGEL+IOFF-IISTA)
#       ifdef _MOLCAS_MPP_
        DBL_MB(MW+IDX-1) = EP
#       else
        GA_Arrays(lg_W)%A(IDX) = EP
#       endif
      end do
    end do

    IOFF = IOFF+NA*NJL
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

iCASE = 7
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
  do ISYA=1,NSYM
    ISYJL = Mul(ISYA,ISYM)
    ISYV = ISYM

    NA = NSSH(ISYA)
    NJL = NIGTJ(ISYJL)
    ! what is start/end in this block?
    IAJGTLSTA = max(IISTA-IOFF,1)
    IAJGTLEND = min(IIEND-IOFF,NA*NJL)

    do IAJGTL=IAJGTLSTA,IAJGTLEND
      IJGTL = (IAJGTL-1)/NA+1
      IA = IAJGTL-NA*(IJGTL-1)
      IJGTLTOT = IJGTL+NIGTJES(ISYJL)
      IJABS = MIGTJ(1,IJGTLTOT)
      ILABS = MIGTJ(2,IJGTLTOT)
      IJ = MIREL(1,IJABS)
      ISYJ = MIREL(2,IJABS)
      IL = MIREL(1,ILABS)
      ISYL = MIREL(2,ILABS)
      do IV=IASTA,IAEND ! these are always all elements
        ! compute integrals (ajvl) and (alvj)
        NV = NVTOT_CHOSYM(Mul(ISYA,ISYJ))
        IAJ = IA-1+NSSH(ISYA)*(IJ-1)
        IVL = IV-1+NASH(ISYV)*(IL-1)
        IOFFAJ = 1+IOBRA(ISYA,ISYJ)+NV*IAJ
        IOFFVL = 1+IOKET(ISYV,ISYL)+NV*IVL
        AJVL = DDOT_(NV,BRABUF(IOFFAJ),1,KETBUF(IOFFVL),1)

        NV = NVTOT_CHOSYM(Mul(ISYA,ISYL))
        IAL = IA-1+NSSH(ISYA)*(IL-1)
        IVJ = IV-1+NASH(ISYV)*(IJ-1)
        IOFFAL = 1+IOBRA(ISYA,ISYL)+NV*IAL
        IOFFVJ = 1+IOKET(ISYV,ISYJ)+NV*IVJ
        ALVJ = DDOT_(NV,BRABUF(IOFFAL),1,KETBUF(IOFFVJ),1)

        ! EM(v,ajl)=((aj,vl)-(al,vj))*SQRT(3/2)
        EM = SQRTA*(AJVL-ALVJ)
        ! write element EM
        IDX = IV+NAS*(IAJGTL+IOFF-IISTA)
#       ifdef _MOLCAS_MPP_
        DBL_MB(MW+IDX-1) = EM
#       else
        GA_Arrays(lg_W)%A(IDX) = EM
#       endif
      end do
    end do

    IOFF = IOFF+NA*NJL
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

call mma_deallocate(BRABUF)
call mma_deallocate(KETBUF)

end subroutine RHSOD_E
