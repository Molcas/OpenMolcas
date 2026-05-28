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

subroutine RHSOD_A_NOSYM(IVEC)

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use SUPERINDEX, only: MTREL, MTUV
use CHOVEC_IO, only: ChoVec_read, ChoVec_Size, NVTOT_ChoSym
use PrintLevel, only: DEBUG
#ifndef _MOLCAS_MPP_
use fake_GA, only: GA_Arrays
#endif
use caspt2_global, only: FIMO, iPrGlb
use caspt2_module, only: nActEl, nAsh, nIsh, nOrb, nSym, NTUV, nTUVES
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) :: iAEnd, iASta, iCase, IDX, IIEnd, IISTA, IJ, IOBRA(8,8), IOFFTJ, IOFFVX, IOKET(8,8), ISYJ, iSym, iSYT, iSYV, &
                     iSYX, IT, ITABS, ITJ, ITTOT, ITVX, ITVXTOT, IV, IVABS, iVX, iX, iXABS, lg_W, mW, nAS, nBra, nFIMOES, nIS, &
                     nKet, nV, nW
real(kind=wp) :: ATVXJ, FTJ, TJVX
real(kind=wp), allocatable :: BRA(:), KET(:)
real(kind=wp), external :: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

if (iPrGlb >= DEBUG) write(u6,*) 'RHS on demand: case A'

!***********************************************************************
! Case A:
!   RHS(tvx,j)=(tj,vx)+FIMO(t,j)*kron(v,x)/NACTEL
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
call CHOVEC_SIZE(1,NBRA,IOBRA)
call CHOVEC_SIZE(2,NKET,IOKET)

call mma_allocate(BRA,NBRA,LABEL='BRA')
call mma_allocate(KET,NKET,LABEL='KET')

call CHOVEC_READ(1,BRA,NBRA)
call CHOVEC_READ(2,KET,NKET)

ICASE = 1
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
NFIMOES = 0
do ISYM=1,NSYM

  NAS = NTUV(ISYM) !NASUP(ISYM,ICASE)
  NIS = NISH(ISYM) !NISUP(ISYM,ICASE)
  NW = NAS*NIS

  if (NW == 0) then
    NFIMOES = NFIMOES+nTri_Elem(NORB(ISYM))
    cycle
  end if

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  do IJ=IISTA,IIEND
    ISYJ = ISYM
    do ITVX=IASTA,IAEND ! these are always all elements
      ITVXTOT = ITVX+NTUVES(ISYM)
      ITABS = MTUV(1,ITVXTOT)
      IVABS = MTUV(2,ITVXTOT)
      IXABS = MTUV(3,ITVXTOT)
      IT = MTREL(1,ITABS)
      ISYT = MTREL(2,ITABS)
      IV = MTREL(1,IVABS)
      ISYV = MTREL(2,IVABS)
      IX = MTREL(1,IXABS)
      ISYX = MTREL(2,IXABS)
      ! compute integrals (tiuv)
      NV = NVTOT_CHOSYM(Mul(ISYT,ISYJ)) ! JSYM=ISYT*ISYI=ISYU*ISYV
      ITJ = IT-1+NASH(ISYT)*(IJ-1)
      IVX = IV-1+NASH(ISYV)*(IX-1)
      IOFFTJ = 1+IOBRA(ISYT,ISYJ)+NV*ITJ
      IOFFVX = 1+IOKET(ISYV,ISYX)+NV*IVX
      TJVX = DDOT_(NV,BRA(IOFFTJ),1,KET(IOFFVX),1)
      ! A(tvx,j) = (tjvx) + FIMO(t,j)*delta(v,x)/NACTEL
      if ((ISYT == ISYJ) .and. (IVABS == IXABS)) then
        ITTOT = IT+NISH(ISYT)
        FTJ = FIMO(NFIMOES+iTri(ITTOT,IJ))
        ATVXJ = TJVX+FTJ/real(max(1,NACTEL),kind=wp)
      else
        ATVXJ = TJVX
      end if
      ! write element A(tvx,j)
      IDX = ITVX+NAS*(IJ-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = ATVXJ
#     else
      GA_Arrays(lg_W)%A(IDX) = ATVXJ
#     endif
    end do
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)

  NFIMOES = NFIMOES+nTri_Elem(NORB(ISYM))

end do
!***********************************************************************

call mma_deallocate(BRA)
call mma_deallocate(KET)

end subroutine RHSOD_A_NOSYM
