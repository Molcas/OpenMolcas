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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

subroutine REO_GASDET(IBLOCK,NBLOCK,ISYM,IREO)
! Create reorder array for determinants : configuration order => Ab order
!
! Jeppe Olsen, November 2001, from GASANA

use lucia_data, only: CONF_REO, IB_CONF_REO, IB_SD_FOR_OPEN, IBCONF_ALL_SYM_FOR_OCCLS, IOCLS, MAXOP, MINOP, MXNSTR, &
                      NCONF_PER_OPEN, NCONF_TOT, NELEC, NGAS, NIRREP, NMXOCCLS, NOBPT, NOCOB, NPDTCNF, NSTSO, NTOOB, PSSIGN
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NBLOCK, IBLOCK(8,NBLOCK), ISYM
integer(kind=iwp), intent(_OUT_) :: IREO(*)
integer(kind=iwp) :: IATP, IBTP, NAEL, NBEL, NEL
integer(kind=iwp), allocatable :: DET_MS(:), DET_OC(:), DET_VC(:), LASTR(:), LBSTR(:), LOCMAX(:), LOCMIN(:), Z(:), ZSCR(:)

!write(u6,*) 'nconf_per_open in reo_gasdet'
!call iwrtma(nconf_per_open,1,4,1,4)

! Specifications of internal space

! Type of alpha and beta strings
IATP = 1
IBTP = 2

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
NEL = NAEL+NBEL

! Info on block structure of space

! Space for alpha and beta strings
call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
! Space for constructing arc weights for configurations
call mma_allocate(ZSCR,(NOCOB+1)*(NEL+1),Label='ZSCR')
call mma_allocate(Z,NOCOB*NEL*2,Label='Z')
call mma_allocate(LOCMIN,NOCOB,Label='LOCMIN')
call mma_allocate(LOCMAX,NOCOB,Label='LOCMAX')
! Occupation and projections of a given determinant
call mma_allocate(DET_OC,NAEL+NBEL,Label='DET_OC')
call mma_allocate(DET_MS,NAEL+NBEL,Label='DET_MS')
call mma_allocate(DET_VC,NOCOB,Label='DET_VC')

!/ Jesper Wisborg Krogh, 2005-06-22
call REO_GASDET_S(IREO,NSTSO(IATP)%A,NSTSO(IBTP)%A,NBLOCK,IBLOCK,NAEL,NBEL,LASTR,LBSTR,NIRREP,NMXOCCLS,NGAS,IOCLS,NTOOB,NOBPT, &
                  IB_CONF_REO,CONF_REO(ISYM)%A,nconf_tot,ib_conf_reo,maxop,nconf_per_open(:,isym),IB_SD_FOR_OPEN,ZSCR,Z,LOCMIN, &
                  LOCMAX,DET_OC,DET_MS,DET_VC,MINOP,IBCONF_ALL_SYM_FOR_OCCLS,PSSIGN,NPDTCNF)

call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(ZSCR)
call mma_deallocate(Z)
call mma_deallocate(LOCMIN)
call mma_deallocate(LOCMAX)
call mma_deallocate(DET_OC)
call mma_deallocate(DET_MS)
call mma_deallocate(DET_VC)

end subroutine REO_GASDET
