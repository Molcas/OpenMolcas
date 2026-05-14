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

subroutine CNFORD_GAS(IOCCLS,NOCCLS,ISYM,ICTSDT,IBLOCK,NBLOCK)
! Generate configurations in ICONF
!
! Generate determinants in configuration order and obtain
! sign array for switching between the two formats.
!
! It is assumed that CSFDIM has been called
!
! Jeppe Olsen Dec. 2001 from CNFORD

use lucia_data, only: CONF_OCC, CONF_REO, IB_CONF_OCC, IB_CONF_REO, IBCONF_ALL_SYM_FOR_OCCLS, MAXOP, MINOP, NCONF_ALL_SYM, &
                      NCONF_PER_OPEN, NCONF_TOT, NGAS, NOBPT, NOCOB
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NOCCLS, IOCCLS(NGAS,NOCCLS), ISYM, NBLOCK, IBLOCK(8,NBLOCK)
integer(kind=iwp), intent(_OUT_) :: ICTSDT(*)
integer(kind=iwp) :: IB_OCCLS, IDOREO, INITIALIZE_CONF_COUNTERS, JOCCLS, NCONF_OCCLS, NCONF_P, NELEC
integer(kind=iwp), allocatable :: LOCMAX(:), LOCMIN(:), Z(:), ZSCR(:)

NELEC = sum(IOCCLS(:,1))

call mma_allocate(ZSCR,(NOCOB+1)*(NELEC+1),Label='ZSCR')
call mma_allocate(Z,NOCOB*NELEC*2,Label='Z')
call mma_allocate(LOCMIN,NOCOB,Label='LOCMIN')
call mma_allocate(LOCMAX,NOCOB,Label='LOCMAX')
! Zero configuration reorder array using NCONF_ALL_SYM
CONF_REO(ISYM)%A(:) = 0

! Generate configurations for all occupation classes

IB_OCCLS = 1
do JOCCLS=1,NOCCLS
  ! Save offset to current occupation class

  ! KIB_OCCLS seems no longer to be used. Therefore the call to
  ! ITOR is commented out to avoid having unitialized arrays
  ! floating around and in case of ITOR there's even written to
  ! this undefined piece of memory.
  !
  ! / Jesper Wisborg Krogh, 2005-06-22

  !_REMOVED call ITOR(WORK(KIB_OCCLS(ISYM)),1,IB_OCCLS,JOCCLS)
  ! Max and min arrays for strings
  call MXMNOC_OCCLS(LOCMIN,LOCMAX,NGAS,NOBPT,IOCCLS(:,JOCCLS),MINOP)
  ! the arcweights
  call CONF_GRAPH(LOCMIN,LOCMAX,NOCOB,NELEC,Z(:),NCONF_P,ZSCR)

  if (JOCCLS == 1) then
    INITIALIZE_CONF_COUNTERS = 1
  else
    INITIALIZE_CONF_COUNTERS = 0
  end if
  IDOREO = 1
  ! Lexical addressing for configurations of this type
  IB_OCCLS = IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS)

  call GEN_CONF_FOR_OCCLS(IOCCLS(:,JOCCLS),IB_OCCLS,INITIALIZE_CONF_COUNTERS,NGAS,ISYM,MINOP,MAXOP,0,NOCOB,NOBPT, &
                          NCONF_PER_OPEN(:,ISYM),NCONF_OCCLS,IB_CONF_REO,IB_CONF_OCC,CONF_OCC(ISYM)%A,IDOREO,Z(:),NCONF_ALL_SYM, &
                          CONF_REO(ISYM)%A,nconf_tot)

  !    GEN_CONF_FOR_OCCLS(IOCCLS,IB_OCCLS,INITIALIZE_CONF_COUNTERS,NGAS,ISYM,MINOP,MAXOP,IONLY_NCONF,NTORB,NOBPT,NCONF_OP, &
  !                       IBCONF_REO,IBCONF_OCC,ICONF,IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM)
  !Error IB_OCCLS = IB_OCCLS+NCONF_ALL_SYM
end do

! Reorder Determinants from configuration order to ab-order

!    REO_GASDET(INBLOCK,NBLOCK,ISYM,IREO,SREO)
call REO_GASDET(IBLOCK,NBLOCK,ISYM,ICTSDT)

!call CNTOST(ICONF,SGNCTS,ICTSDT,IWORK,IDCNF(1),IREFSM,IREFML,IREFG,XNDXCI,NORB,NEL,ORBSYM,ICNFOK,IGENSG,ISGNA,ISGNB,IPRCSF)

call mma_deallocate(ZSCR)
call mma_deallocate(Z)
call mma_deallocate(LOCMIN)
call mma_deallocate(LOCMAX)

end subroutine CNFORD_GAS
