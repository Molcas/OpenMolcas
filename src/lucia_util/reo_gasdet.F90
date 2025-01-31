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

use stdalloc, only: mma_allocate, mma_deallocate
use GLBBAS, only: DFTP, CONF_REO
use strbas, only: NSTSO, IOCLS
use lucia_data, only: IBCONF_ALL_SYM_FOR_OCCLS, IB_CONF_REO, IB_SD_FOR_OPEN, MAXOP, MINOP, NCONF_PER_OPEN, NCONF_TOT, NPDTCNF
use lucia_data, only: NGAS, NMXOCCLS
use lucia_data, only: IPRDIA
use lucia_data, only: PSSIGN
use lucia_data, only: MXNSTR, IBSPGPFTP, NELFSPGP
use lucia_data, only: NOCOB, NOBPT, NTOOB
use lucia_data, only: NOCTYP
use lucia_data, only: NELEC
use lucia_data, only: MXPNGAS
use csm_data, only: NSMST

implicit none
! =====
! Input
! =====
integer NBLOCK, ISYM
integer IBLOCK(8,NBLOCK)
! Output
integer IREO(*)
integer, allocatable :: LASTR(:), LBSTR(:)
integer, allocatable :: ZSCR(:), Z(:)
integer, allocatable :: LOCMIN(:), LOCMAX(:)
integer, allocatable :: DET_OC(:), DET_MS(:), DET_VC(:)
integer NTEST, IATP, IBTP, NAEL, NBEL, NEL, NOCTPA, NOCTPB, IOCTPA, IOCTPB

!write(6,*) 'nconf_per_open in reo_gasdet'
!call iwrtma(nconf_per_open,1,4,1,4)

! Specifications of internal space

NTEST = 0
NTEST = max(NTEST,IPRDIA)
! Type of alpha and beta strings
IATP = 1
IBTP = 2

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
NEL = NAEL+NBEL

NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)

IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)

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
call REO_GASDET_S(IREO,NSTSO(IATP)%I,NSTSO(IBTP)%I,NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,NAEL,NBEL,LASTR,LBSTR,NSMST, &
                  NELFSPGP,NMXOCCLS,NGAS,IOCLS,NTOOB,NOBPT,DFTP,IB_CONF_REO,conf_reo(isym)%I,nconf_tot,ib_conf_reo,maxop, &
                  nconf_per_open(1,isym),IB_SD_FOR_OPEN,ZSCR,Z,LOCMIN,LOCMAX,DET_OC,DET_MS,DET_VC,MINOP,IBCONF_ALL_SYM_FOR_OCCLS, &
                  PSSIGN,NPDTCNF)

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
