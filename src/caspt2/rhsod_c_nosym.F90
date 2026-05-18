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

subroutine RHSOD_C_NOSYM(IVEC)

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp
use Constants, only: Zero
use SUPERINDEX, only: MTUV, MTREL, KTUV
use CHOVEC_IO, only: NVTOT_CHOSYM, CHOVEC_SIZE, CHOVEC_READ
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG
use caspt2_global, only: FIMO
use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
use fake_GA, only: GA_Arrays
#endif
use caspt2_module, only: nActEl, nAshT, nSym, nASup, nISup, NTUVES, nSsh, nAsh, nIsh, NAES, nOrb

implicit none
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) IOBRA(8,8), IOKET(8,8)
real(kind=wp), allocatable :: BRA(:), KET(:)
real(kind=wp) AddOne, ATVX, FAT, SUMU
real(kind=wp), external :: DDot_
integer(kind=iwp) IA, IAEND, IASTA, IAT, IATOT, iCASE, IDX, IIEND, IISTA, IOFFAT, IOFFVX, ISYA, iSym, ISYT, ISYV, ISYX, IT, ITABS, &
                  ITTOT, ITVV, ITVX, ITVXTOT, IUABS, IUUT, IV, IVABS, IVX, IX, IXABS, lg_W, mW, NAS, NBRA, NFIMOES, NIS, NKET, NV, &
                  NW
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

if (iPrGlb >= DEBUG) write(6,*) 'RHS on demand: case C'

!***********************************************************************
! Case C:
!   RHS(tvx,a)=(at,vx)+(FIMO(a,t)-Sum_u(au,ut))*delta(v,x)/NACTEL
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
call CHOVEC_SIZE(3,NBRA,IOBRA)
call CHOVEC_SIZE(2,NKET,IOKET)

call mma_allocate(BRA,NBRA,LABEL='BRA')
call mma_allocate(KET,NKET,LABEL='KET')

call CHOVEC_READ(3,BRA,NBRA)
call CHOVEC_READ(2,KET,NKET)

ICASE = 4
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
NFIMOES = 0
do ISYM=1,NSYM

  NAS = NASUP(ISYM,ICASE) !NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE) !NISUP(ISYM,ICASE)
  NW = NAS*NIS

  if (NW == 0) then
    NFIMOES = NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2
    cycle
  end if

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  do IA=IISTA,IIEND
    ISYA = ISYM
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
      ! compute integrals (at,vx)
      NV = NVTOT_CHOSYM(Mul(ISYA,ISYT)) ! JSYM=ISYT*ISYI=ISYU*ISYV
      IAT = IA-1+NSSH(ISYA)*(IT-1)
      IVX = IV-1+NASH(ISYV)*(IX-1)
      IOFFAT = 1+IOBRA(ISYA,ISYT)+NV*IAT
      IOFFVX = 1+IOKET(ISYV,ISYX)+NV*IVX
      ATVX = DDOT_(NV,BRA(IOFFAT),1,KET(IOFFVX),1)

      ! W(tvx,a) = (at,vx) + (FIMO(a,t)-Sum_u(au,ut))*delta(v,x)/NACTEL
      ! write element W(tvx,j), only the (at,vx) part
      IDX = ITVX+NAS*(IA-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = ATVX
#     else
      GA_Arrays(lg_W)%A(IDX) = ATVX
#     endif
    end do
    ! now, add in the part with corrections to the integrals
    IATOT = IA+NISH(ISYM)+NASH(ISYM)
    do IT=1,NASH(ISYM)
      ITTOT = IT+NISH(ISYM)
      FAT = FIMO(NFIMOES+(IATOT*(IATOT-1))/2+ITTOT)
      SUMU = Zero
      ITABS = NAES(ISYM)+IT
      do IUABS=1,NASHT
        IUUT = KTUV(IUABS,IUABS,ITABS)-NTUVES(ISYM)
        IDX = IUUT+NAS*(IA-IISTA)
#       ifdef _MOLCAS_MPP_
        SUMU = SUMU+DBL_MB(MW+IDX-1)
#       else
        SUMU = SUMU+GA_Arrays(lg_w)%A(IDX)
#       endif
      end do
      ADDONE = (FAT-SUMU)/dble(max(1,NACTEL))
      do IVABS=1,NASHT
        ITVV = KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
        IDX = ITVV+NAS*(IA-IISTA)
#       ifdef _MOLCAS_MPP_
        DBL_MB(MW+IDX-1) = DBL_MB(MW+IDX-1)+ADDONE
#       else
        GA_Arrays(lg_w)%A(IDX) = GA_Arrays(lg_w)%A(IDX)+ADDONE
#       endif
      end do
    end do
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)

  NFIMOES = NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

end do
!***********************************************************************

call mma_deallocate(BRA)
call mma_deallocate(KET)

end subroutine RHSOD_C_NOSYM
