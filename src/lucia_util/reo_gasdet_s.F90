!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine REO_GASDET_S(IREO,NSSOA,NSSOB,NBLOCK,IBLOCK,NAEL,NBEL,IASTR,IBSTR,NSMST,NOCCLS,NGAS,IOCCLS,NORB,NOBPT,IB_CONF_OPEN, &
                        iconf_reo,nconf_tot,ib_conf_reo,maxop,nconf_per_open,IB_SD_FOR_OPEN,IZSCR,IZ,IOCMIN,IOCMAX,IDET_OC, &
                        IDET_MS,IDET_VC,MINOP,IBCONF_ALL_SYM_FOR_OCCLS,PSSIGN,NPDTCNF)
! SUBROUTINE REO_GASDET_S --> 44
!
! Reorder determinants in GAS space from det to configuration order
!
! IBCONF_ALL_SYM_FOR_OCCLS : Offset to start of configurations of given occls in list containing all symmetries
! Z_PTDT(IOPEN+1)%A gives Z array for prototype dets with IOPEN
! REO_PTDT(IOPEN+1)%A gives the corresponding reorder array open orbitals

use lucia_data, only: REO_PTDT, Z_PTDT
use Constants, only: One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: IREO(*), IZSCR(*), IZ(*), IOCMIN(*), IOCMAX(*), IDET_OC(*), IDET_MS(*), IDET_VC(*)
integer(kind=iwp), intent(in) :: NSMST, NSSOA(NSMST,*), NSSOB(NSMST,*), NBLOCK, IBLOCK(8,NBLOCK), NAEL, NBEL, NOCCLS, NGAS, &
                                 IOCCLS(NGAS,NOCCLS), NORB, NOBPT(*), IB_CONF_OPEN(*), nconf_tot, iconf_reo(nconf_tot), maxop, &
                                 ib_conf_reo(maxop+1), nconf_per_open(maxop+1), IB_SD_FOR_OPEN(*), MINOP, &
                                 IBCONF_ALL_SYM_FOR_OCCLS(NOCCLS), NPDTCNF(*)
integer(kind=iwp), intent(inout) :: IASTR(*), IBSTR(*)
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp) :: IA, IADR_SD_CONF_ORDER, IAGRP, IASM, IATP, IB, IB_OCCLS, IBCNF_OUT, IBGRP, IBSM, IBTP, icnf_out, IDET, &
                     IDUM(1), IOC, IPTDT, IRESTR, ISIGN_2003, ISGN, JBLOCK, MINIA, NASTR1, NBSTR1, nconf_op, NCONF_P, NDOUBLE, &
                     NEL, NIA, NIB, NOCOB, NOPEN, NOPEN_AL, NPTDT
integer(kind=iwp), external :: ilex_for_conf_new, IZNUM_PTDT, NOP_FOR_CONF

IAGRP = 1
IBGRP = 2

NEL = NAEL+NBEL

IDET = 0
do JBLOCK=1,NBLOCK
  IATP = IBLOCK(1,JBLOCK)
  IBTP = IBLOCK(2,JBLOCK)
  IASM = IBLOCK(3,JBLOCK)
  IBSM = IBLOCK(4,JBLOCK)
  !write(u6,*) ' REO_GASDET, IATP, IBTP = ',IATP,IBTP
  ! Occupation class of this combination of string
  call IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
  !    IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
  ! Arcweights for this occupation class
  call MXMNOC_OCCLS(IOCMIN,IOCMAX,NGAS,NOBPT,IOCCLS(:,IOC),MINOP)
  !    MXMNOC_OCCLS(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP)
  ! the arcweights
  call CONF_GRAPH(IOCMIN,IOCMAX,NORB,NEL,IZ,NCONF_P,IZSCR)
  !    CONF_GRAPH(IOCC_MIN,IOCC_MAX,NORB,NEL,IARCW,NCONF,ISCR)
  ! Obtain alpha strings of sym IASM and type IATP
  IDUM(1) = 0
  call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,NORB,0,IDUM,IDUM)
  ! Obtain Beta strings of sym IBSM and type IBTP
  IDUM(1) = 0
  call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,NORB,0,IDUM,IDUM)
  ! Occupation class corresponding to this combination
  ! The following call should presumably use 'IOC' rather than 'IOCNUM'
  ! The variable name IOCNUM seems to be used nowhere... PAM 2009
  !call IAIB_TO_OCCLS(1,IATP,2,IBTP,IOCNUM)
  call IAIB_TO_OCCLS(1,IATP,2,IBTP,IOC)
  !    IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
  ! Offset to this occupation class in occupation class ordered cnf list
  IB_OCCLS = IBCONF_ALL_SYM_FOR_OCCLS(IOC)
  ! Info for this occupation class:
  IRESTR = 0
  if ((PSSIGN == One) .and. (IASM == IBSM) .and. (IATP == IBTP)) IRESTR = 1

  NIA = NSSOA(IASM,IATP)
  NIB = NSSOB(IBSM,IBTP)

  do IB=1,NIB
    if (IRESTR == 1) then
      MINIA = IB
    else
      MINIA = 1
    end if
    do IA=MINIA,NIA
      IDET = IDET+1
      !PAM06 call ABSTR_TO_ORDSTR(IASTR(1,IA),IBSTR(1,IB),NAEL,NBEL,
      call ABSTR_TO_ORDSTR(IASTR(1+NAEL*(IA-1)),IBSTR(1+NBEL*(IB-1)),NAEL,NBEL,IDET_OC,IDET_MS,ISGN)
      !    ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,IDET_OC,IDET_SP,ISGN)
      ! Number of open orbitals in this configuration
      NOPEN = NOP_FOR_CONF(IDET_OC,NEL)
      !       NOP_FOR_CONF(ICONF,NEL)
      NDOUBLE = (NEL-NOPEN)/2
      NOCOB = NOPEN+NDOUBLE
      NOPEN_AL = NAEL-NDOUBLE
      !write(u6,*) ' NOPEN, NOPEN_AL = ',NOPEN,NOPEN_AL
      !ERROR NPTDT = IBINOM(NOPEN,NOPEN_AL)
      NPTDT = NPDTCNF(NOPEN+1)
      ! Packed form of this configuration
      call REFORM_CONF_OCC(IDET_OC,IDET_VC,NEL,NOCOB)
      ! Address of this configuration
      ! Offset to configurations with this number of open orbitals in
      ! reordered cnf list
      !write(u6,*) 'iconf_reo_new array:'
      !call iwrtma(iconf_reo_new,1,nconf_tot,1,nconf_tot)
      !.. Giovanni and Dongxia comment off the following line
      !ICNF_OUT = ILEX_FOR_CONF(IDET_VC,NOCOB,NORB,NEL,IZ,1,ICONF_REO(IB_OCCLS))
      !!          ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
      !.. end
      !write(u6,*) 'ib_conf_reo at line 2401, and maxop',maxop
      !call iwrtma(ib_conf_reo,1,maxop+1,1,maxop+1)
      !call iwrtma(nconf_per_open,1,maxop+1,1,maxop+1)
      !write(u6,*) 'before calling ilex_for_conf_new, in reogas_det_s'
      !write(u6,*) 'nopen =',nopen
      !write(u6,*) 'and nconf_per_open(nopen+1) =',nconf_per_open(nopen+1)
      !write(u6,*) 'check iconf_reo array'
      !call iwrtma(iconf_reo,1,nconf_tot,1,nconf_tot)
      nconf_op = nconf_per_open(nopen+1)
      !call iwrtma(nconf_per_open,1,maxop+1,1,maxop+1)
      icnf_out = ilex_for_conf_new(idet_vc,nocob,norb,nel,iz,1,iconf_reo(ib_conf_reo(nopen+1)),nconf_op,ib_occls)+ &
                 ib_conf_reo(nopen+1)-1
      !write(u6,*) ' number of configuration in output list',ICNF_OUT
      ! Spinprojections of open orbitals
      call EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_VC,NEL)
      !    EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_OPEN_MS,NEL)

      ISIGN_2003 = 1
      if (abs(PSSIGN) == One) then
        ! If combinations are used, then the prototype determinants
        ! are defined so the first open spin-orbital is having alpha spin.
        ! In ab order, the included determinant is defined, by having
        ! alpha-spin in the first singly occupied orbital. These definitions
        ! may differ, so ensure that the included det obeys prototype constraint
        ! Address of this spinprojection pattern
        if (IDET_VC(1) < 0) then
          IDET_VC(1:NOPEN) = -IDET_VC(1:NOPEN)
          if (PSSIGN == -One) ISIGN_2003 = -1
          ! Update sign AB => ordered list
          !PAM06 call ABSTR_TO_ORDSTR(IBSTR(1,IB),IASTR(1,IA),NBEL,NAEL,
          call ABSTR_TO_ORDSTR(IBSTR(1+NBEL*(IB-1)),IASTR(1+NAEL*(IA-1)),NBEL,NAEL,IDET_OC,IDET_MS,ISGN)
        end if
      end if
      IPTDT = IZNUM_PTDT(IDET_VC,NOPEN,NOPEN_AL,Z_PTDT(NOPEN+1)%A,REO_PTDT(NOPEN+1)%A,1)
      !       IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
      !write(u6,*) ' Number of det in list of PTDT ', IPTDT
      !write(u6,*) ' IB_SD_FOR_OPEN(NOPEN+1) = ',IB_SD_FOR_OPEN(NOPEN+1)
      !write(u6,*) ' ICNF_OUT, NPTDT ',ICNF_OUT, NPTDT
      IBCNF_OUT = IB_CONF_OPEN(NOPEN+1)
      !write(u6,*) ' IBCNF_OUT = ',IBCNF_OUT
      IADR_SD_CONF_ORDER = IB_SD_FOR_OPEN(NOPEN+1)-1+(ICNF_OUT-IBCNF_OUT)*NPTDT+IPTDT
      if (IADR_SD_CONF_ORDER <= 0) then
        write(u6,*) ' Problemo, IADR_SD_CONF_ORDER < 0'
        write(u6,*) ' IADR_SD_CONF_ORDER = ',IADR_SD_CONF_ORDER
        call XFLUSH(u6)
      end if
      !write(u6,*) ' IADR_SD_CONF_ORDER, ISGN, IDET = ',IADR_SD_CONF_ORDER,ISGN,IDET
      IREO(IADR_SD_CONF_ORDER) = ISGN*IDET*ISIGN_2003

    end do
    ! End of loop over alpha strings
  end do
  ! End of loop over beta strings
end do
! End of loop over blocks

#ifdef _DEBUGPRINT_
write(u6,*) ' Reorder array, CONF order => string order'
write(u6,*) ' ========================================='
call IWRTMA(IREO,1,IDET,1,IDET)
#endif

end subroutine REO_GASDET_S
