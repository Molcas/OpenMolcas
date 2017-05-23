************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE REO_GASDET_S(   IREO,  NSSOA,  NSSOB, NOCTPA, NOCTPB,
     &                        MXPNGAS, IOCTPA, IOCTPB, NBLOCK, IBLOCK,
     &                           NAEL,   NBEL,  IASTR,  IBSTR,
     &                          NSMST,NELFSPGP,NOCCLS,  NGAS,IOCCLS,
     &                           NORB,NOBPT,DFTP,IB_CONF_OPEN,iconf_reo,
*
     &                        nconf_tot,
     &                        ib_conf_reo,
     &                          maxop,
     &                        nconf_per_open,
     &                        IB_SD_FOR_OPEN,IZSCR,IZ,IOCMIN,IOCMAX,
*
     &                        IDET_OC,IDET_MS,IDET_VC,  iWORK,KZ_PTDT,
     &                        KREO_PTDT,
     &                          MINOP,
     &                        IBCONF_ALL_SYM_FOR_OCCLS,
     &                         PSSIGN,
*
     &                        NPDTCNF)
*
* SUBROUTINE REO_GASDET_S --> 44
*
*
* Reorder determinants in GAS space from det to configuration order
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
      INTEGER NOBPT(*)
      INTEGER DFTP(*)
      INTEGER IB_CONF_OPEN(*),iconf_reo(nconf_tot)
      integer ib_conf_reo(maxop+1),nconf_per_open(maxop+1)
      INTEGER IB_SD_FOR_OPEN(*)
      INTEGER NPDTCNF(*)
*. Offset to start of configurations of given occls in list containing all symmetries
      INTEGER IBCONF_ALL_SYM_FOR_OCCLS(NOCCLS)
*. iWORK(KZ_PTDT(IOPEN+1) gives Z  array for prototype dets with IOPEN
*. iWORK(KREO_PTDT(IOPEN+1) gives the corresponding reorder array
*. open orbitals
      INTEGER KZ_PTDT(*), KREO_PTDT(*)
*. The work array containing used for iWORK(KZ_PTDET()),iWORK(KREO_PTDT())
      DIMENSION iWORK(*)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
*. Scratch space
*PAM06 Note: NAEL and NBEL could legally be =0!
*PAM06      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IASTR(*),IBSTR(*)
      INTEGER IZSCR(*),IZ(*),IOCMIN(*),IOCMAX(*)
      INTEGER IDET_OC(*), IDET_MS(*) , IDET_VC(*)
*. Output
      INTEGER IREO(*)
C     DIMENSION SREO(*)
*
      NTEST = 0
*
      IAGRP = 1
      IBGRP = 2
*
      NEL = NAEL + NBEL
*
      IDET = 0
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
C?      WRITE(6,*) ' REO_GASDET, IATP, IBTP = ', IATP, IBTP
*. Occupation class of this combination of string
        CALL IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
C            IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
*. Arcweights for this occupation class
        CALL MXMNOC_OCCLS(IOCMIN,IOCMAX,    NGAS,   NOBPT,IOCCLS(1,IOC),
     &                       MINOP,   NTEST)
C     MXMNOC_OCCLS(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,NTESTG)
*. the arcweights
         CALL CONF_GRAPH(  IOCMIN,  IOCMAX,    NORB,     NEL,      IZ,
     &                    NCONF_P,   IZSCR)
C             CONF_GRAPH(IOCC_MIN,IOCC_MAX,NORB,NEL,IARCW,NCONF,ISCR)
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(      1,   IATP,   IASM,   NAEL, NASTR1,
     &                           IASTR,   NORB,      0,   IDUM,   IDUM)
*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(      2,   IBTP,   IBSM,   NBEL, NBSTR1,
     &                           IBSTR,   NORB,      0,   IDUM,   IDUM)
*. Occupation class corresponding to this combination
C            IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
* The following call should presumably use 'IOC' rather than 'IOCNUM'
* The variable name IOCNUM seems to be used nowhere... PAM 2009
*        CALL IAIB_TO_OCCLS(1,IATP,2,IBTP,IOCNUM)
        CALL IAIB_TO_OCCLS(1,IATP,2,IBTP,IOC)
*. Offset to this occupation class in occupation class ordered cnf list
        IB_OCCLS = IBCONF_ALL_SYM_FOR_OCCLS(IOC)
*. Info for this occupation class :
        IRESTR = 0
        IF(PSSIGN.EQ.1.0D0.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
         IRESTR = 1
        END IF
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        DO  IB = 1,NIB
          IF(IRESTR.EQ.1) THEN
            MINIA = IB
          ELSE
            MINIA = 1
          END IF
          DO  IA = MINIA,NIA
            IDET = IDET + 1
C                ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,IDET_OC,IDET_SP,ISIGN)
*PAM06       CALL ABSTR_TO_ORDSTR(IASTR(1,IA),IBSTR(1,IB),NAEL,NBEL,
            CALL ABSTR_TO_ORDSTR(IASTR(1+NAEL*(IA-1)),
     &                           IBSTR(1+NBEL*(IB-1)),
     &                             NAEL,  NBEL,IDET_OC,IDET_MS,ISIGN)
*. Number of open orbitals in this configuration
            NOPEN = NOP_FOR_CONF(IDET_OC,NEL)
C                   NOP_FOR_CONF(ICONF,NEL)
            NDOUBLE = (NEL-NOPEN)/2
            NOCOB = NOPEN + NDOUBLE
            NOPEN_AL = NAEL - NDOUBLE
C?          WRITE(6,*) ' NOPEN, NOPEN_AL = ', NOPEN,NOPEN_AL
CERRROR     NPTDT = IBION_LUCIA(NOPEN,NOPEN_AL)
            NPTDT = NPDTCNF(NOPEN+1)
*. Packed form of this configuration
C                REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY)
            CALL REFORM_CONF_OCC(IDET_OC,IDET_VC,NEL,NOCOB,1)
*. Address of this configuration
*. Offset to configurations with this number of open orbitals in
*. reordered cnf list
C                      ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
c           write(6,*)'iconf_reo_new array:'
c           call iwrtma(iconf_reo_new,1,nconf_tot,1,nconf_tot)
c.... Giovanni and Dongxia comment off the following 2 lines
c           ICNF_OUT = ILEX_FOR_CONF(IDET_VC,NOCOB,NORB,NEL,IZ,1,
c    &                 ICONF_REO(IB_OCCLS))
c..... end
c           write(6,*)'ib_conf_reo at line 2401, and maxop',maxop
c           call iwrtma(ib_conf_reo,1,maxop+1,1,maxop+1)
c           call iwrtma(nconf_per_open,1,maxop+1,1,maxop+1)
c          write(6,*)'before calling ilex_for_conf_new, in reogas_det_s'
c          write(6,*)'nopen =',nopen
c          write(6,*)'and nconf_per_open(nopen+1) =',
c    &          nconf_per_open(nopen+1)
c          write(6,*)'check iconf_reo array'
c          call iwrtma(iconf_reo,1,nconf_tot,1,nconf_tot)
            nconf_op = nconf_per_open(nopen+1)
c          call iwrtma(nconf_per_open,1,maxop+1,1,maxop+1)
            icnf_out=ilex_for_conf_new(idet_vc,nocob,norb,nel,iz,1,
     &          iconf_reo(ib_conf_reo(nopen+1)),nconf_op,ib_occls)
     &         +ib_conf_reo(nopen+1)-1
C?          WRITE(6,*) ' number of configuration in output list',
C?   &      ICNF_OUT
*. Spinprojections of open orbitals
            CALL EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_VC,NEL)
C                EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_OPEN_MS,NEL)

            ISIGN_2003 = 1
            IF(ABS(PSSIGN).EQ.1.0D0) THEN
*. If combinations are used, then the prototype determinants
*. are defined so the first open spin-orbital is having alpha spin.
*. In ab order, the included determinant is defined, by having
*. alpha-spin in the first singly occupied orbital. These definitions
*. may differ, so ensure that the included det obeys prototype constraint
*. Address of this spinprojection pattern
              IF(IDET_VC(1).LT.0) THEN
                DO I = 1, NOPEN
                  IDET_VC(I) = -1*IDET_VC(I)
                END DO
                IF(PSSIGN.EQ.-1.0D0) ISIGN_2003 = -1
*. Update sign AB => ordered list
*PAM06          CALL ABSTR_TO_ORDSTR(IBSTR(1,IB),IASTR(1,IA),NBEL,NAEL,
            CALL ABSTR_TO_ORDSTR(IBSTR(1+NBEL*(IB-1)),
     &                           IASTR(1+NAEL*(IA-1)),
     &                             NBEL,  NAEL,IDET_OC,IDET_MS,ISIGN)
              END IF
           END IF
C  IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
            IPTDT = IZNUM_PTDT(IDET_VC,NOPEN,NOPEN_AL,
     &              iWORK(KZ_PTDT(NOPEN+1)),iWORK(KREO_PTDT(NOPEN+1)),
     &              1)
C?          WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
C?          WRITE(6,*) ' IB_SD_FOR_OPEN(NOPEN+1) = ',
C?   &                   IB_SD_FOR_OPEN(NOPEN+1)
C?          WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
            IBCNF_OUT = IB_CONF_OPEN(NOPEN+1)
C?          WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
            IADR_SD_CONF_ORDER = IB_SD_FOR_OPEN(NOPEN+1) - 1
     &                         + (ICNF_OUT-IBCNF_OUT)*NPTDT + IPTDT
            IF(IADR_SD_CONF_ORDER.LE.0) THEN
              WRITE(6,*) ' Problemo, IADR_SD_CONF_ORDER < 0 '
              WRITE(6,*) ' IADR_SD_CONF_ORDER = ', IADR_SD_CONF_ORDER
              CALL XFLUSH(6)
            END IF
C?          WRITE(6,*) ' IADR_SD_CONF_ORDER, ISIGN, IDET = ',
C?   &                   IADR_SD_CONF_ORDER, ISIGN, IDET
            IREO(IADR_SD_CONF_ORDER) = ISIGN*IDET*ISIGN_2003
*
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
        END DO
*       ^ End of loop over blocks
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reorder array, CONF order => string order '
        WRITE(6,*) ' ========================================== '
        CALL IWRTMA(IREO,1,IDET,1,IDET)
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NOCTPA)
        CALL Unused_integer(NOCTPB)
        CALL Unused_integer(IOCTPA)
        CALL Unused_integer(IOCTPB)
        CALL Unused_integer_array(NELFSPGP)
        CALL Unused_integer_array(DFTP)
      END IF
      END
*
