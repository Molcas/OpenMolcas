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
      SUBROUTINE DENSI2_LUCIA(    I12,   RHO1,   RHO2,  RHO2S,  RHO2A,
     &                              L,      R,    LUL,    LUR,  EXPS2,
     &                        IDOSRHO1, SRHO1, IPACK)
*
* Density matrices between L and R
*
* I12 = 1 => only one-body density
* I12 = 2 => one- and two-body density matrices
*
* Jeppe Olsen,      Oct 94
* GAS modifications Aug 95
* Two body density added, '96
*
* Table-Block driven, June 97
* Spin density added, Jan. 99
*
* Jesper Wisborg Krogh
* Allowing to symmetry pack on the fly, Sept. 2003
*
* Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
* ijkl = ij*(ij-1)/2+kl, ij.ge.kl
*
* Two-body symmetric density stored in rho2s
* Two-body anti-symmetric density stored in rho2a
*
* If the two-body density matrix is calculated, then also the
* expectation value of the spin is evaluated.
* The latter is realized as
* S**2
*      = S+S- + Sz(Sz-1)
*      = -Sum(ij) a+i alpha a+j beta a i beta a j alpha + Nalpha +
*        1/2(N alpha - N beta))(1/2(N alpha - Nbeta) - 1)
* If IDOSRHO1 = 1, spin density is also calculated
      IMPLICIT REAL*8(A-H,O-Z)
c      REAL*8 INPRDD
*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "cicisp.fh"
#include "strbas.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "WrkSpc.fh"
#include "crun.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "cprnt.fh"
#include "spinfo_lucia.fh"
#include "glbbas.fh"
*
      LOGICAL IPACK
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
#include "lucinp.fh"
#include "clunit.fh"
*. Scratch for string information
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      INTEGER SXSTSM(1)
*. Specific input
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*),RHO2S(*),RHO2A(*),SRHO1(*)
*. Before I forget it :
      IDUM = 0
*     CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
      ZERO = 0.0D0
      CALL SETVEC(RHO1,ZERO ,NACOB ** 2 )
      IF(I12.EQ.2) THEN
         IF(IPACK) THEN
* If IPACK .EQ. .TRUE. then
C     Number of elements in symmetric and antisymmetric 2-body
C     density matrices are given in Nijkl.
            NIJ   = (NACOB*(NACOB+1))/2
            NIJKL = (NIJ*(NIJ+1))/2
            CALL SETVEC(RHO2S,ZERO,NIJKL)
            CALL SETVEC(RHO2A,ZERO,NIJKL)
         ELSE
            CALL SETVEC(RHO2,ZERO ,NACOB ** 2 *(NACOB**2+1)/2)
         END IF
      END IF
*
      IF(IDOSRHO1.EQ.1) THEN
        CALL SETVEC(SRHO1,ZERO,NACOB ** 2)
      END IF
*
C?     WRITE(6,*) ' ISSPC ICSPC in DENSI2 ',ISSPC,ICSPC
*
* Info for this internal space
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
      JATP = 1
      JBTP = 2
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offsets for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      ILSM = ISSM
      IRSM = ICSM

* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL GETMEM('KSTSTS','ALLO','INTE',KSTSTS,NSMST ** 2)
      CALL GETMEM('KSTSTD','ALLO','INTE',KSTSTD,NSMST ** 2)
      CALL STSTSM(IWORK(KSTSTS),IWORK(KSTSTD),NSMST)
*. connection matrices for supergroups
      CALL GETMEM('CONSPA','ALLO','INTE',KCONSPA,NOCTPA**2)
      CALL GETMEM('CONSPB','ALLO','INTE',KCONSPB,NOCTPB**2)
      CALL SPGRPCON(   IOCTPA,   NOCTPA,     NGAS,  MXPNGAS, NELFSPGP,
     &              iWORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(   IOCTPB,   NOCTPB,     NGAS,  MXPNGAS, NELFSPGP,
     &              iWORK(KCONSPB),IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(IWORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(IWORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(IWORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(IWORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(IWORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(IWORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MAXA = MAX(MAXA,MAXA0)
      MAXB = MAX(MAXB,MAXB0)
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRDEN.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
*. replace with MXINKA !!!
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
C?    WRITE(6,*) ' DENSI2 : MAXI MAXK ', MAXI,MAXK
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
         DO IOBSM = 1, NSMOB
            MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
         END DO
      END DO
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(IPRDEN.GE.2) write(6,*) ' DENSI2 : MXSB MXTSOB MXSOOB ',
     &       MXSB,MXTSOB,MXSOOB
c      IF(ISIMSYM.NE.1) THEN
        LSCR1 = MXSOOB
c      ELSE
c        LSCR1 = MXSOOB_AS
c      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
* JESPER: Should reduce I/O
      IF (ENVIRO(1:6).EQ.'RASSCF') THEN
        LSCR1 = MAX(INT(XISPSM(IREFSM, 1)),MXSOOB)
        IF(PSSIGN.NE.0.0D0) LSCR1 = INT(2.0D0*XISPSM(IREFSM,1))
      ENDIF
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      CALL GETMEM('INSCR ','ALLO','REAL',KINSCR,INTSCR)
*
*. Arrays giving allowed type combinations '
      CALL GETMEM('SIOIO ','ALLO','INTE',KSIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('CIOIO ','ALLO','INTE',KCIOIO,NOCTPA*NOCTPB)
*
      CALL IAIBCM(ISSPC,iWORK(KSIOIO))
      CALL IAIBCM(ISSPC,iWORK(KCIOIO))
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESCPH(iWORK(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                  NSMST,NSTFSMSPGP,MXPNSMST,   NSMOB, MXPNGAS,
     &                   NGAS,   NOBPTS,   IPRCIX,     MAXK, NELFSPGP,
     &                   MXCJ,   MXCIJA,   MXCIJB,  MXCIJAB,   MXSXBL,
     &               MXADKBLK,   IPHGAS, NHLFSPGP,     MNHL,  IADVICE,
*
     &              MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
      IF(IPRDEN.GE.2) THEN
        WRITE(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
      LSCR12 = MAX(LSCR1,2*LSCR2)
      IF (ENVIRO(1:6) .EQ. 'RASSCF') THEN
         LSCR12 = MAX(LSCR1,LSCR2)
      END IF
*. It is assumed that the third block already has been allocated, so
      KC2 = KVEC3
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
      IF (ENVIRO(1:6) .EQ. 'RASSCF') THEN
         KCSCR = KC2 + LSCR12
      END IF
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL GETMEM('I1    ','ALLO','INTE',KI1,  LSCR3       )
      CALL GETMEM('I2    ','ALLO','INTE',KI2,  LSCR3       )
      CALL GETMEM('I3    ','ALLO','INTE',KI3,  LSCR3       )
      CALL GETMEM('I4    ','ALLO','INTE',KI4,  LSCR3       )
      CALL GETMEM('XI1S  ','ALLO','REAL',KXI1S,LSCR3       )
      CALL GETMEM('XI2S  ','ALLO','REAL',KXI2S,LSCR3       )
      CALL GETMEM('XI3S  ','ALLO','REAL',KXI3S,LSCR3       )
      CALL GETMEM('XI4S  ','ALLO','REAL',KXI4S,LSCR3       )
*. Arrays giving block type
      CALL GETMEM('SBLTP ','ALLO','INTE',KSBLTP,NSMST)
      CALL GETMEM('CBLTP ','ALLO','INTE',KCBLTP,NSMST)
*. Arrays for additional symmetry operation
c      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
c        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
c        CALL SIGVST(WORK(KSVST),NSMST)
c      ELSE
         KSVST = 1
c      END IF
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,iWORK(KSBLTP),iWORK(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,iWORK(KCBLTP),iWORK(KSVST))
*.0 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
* scratch space containing active one body
      CALL GETMEM('RHO1S ','ALLO','REAL',KRHO1S,NACOB ** 2)
*. For natural orbitals
      CALL GETMEM('RHO1P ','ALLO','REAL',KRHO1P,NACOB*(NACOB+1)/2)
      CALL GETMEM('XNATO ','ALLO','REAL',KXNATO,NACOB **2)
*. Natural orbitals in symmetry blocks
      CALL GETMEM('RHO1S ','ALLO','REAL',KRHO1SM,NACOB ** 2)
      CALL GETMEM('RHO1S ','ALLO','REAL',KXNATSM,NACOB ** 2)
      CALL GETMEM('RHO1S ','ALLO','REAL',KOCCSM,NACOB )
*
*. Space for one block of string occupations and two arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      CALL GETMEM('KLZSCR','ALLO','INTE',KLZSCR,LZSCR)
      DO K12 = 1, 1
        CALL GETMEM('KLOCS ','ALLO','INTE',KLOCSTR(K12),MAX_STR_OC_BLK)
      END DO
      DO I1234 = 1, 2
        CALL GETMEM('KLREO ','ALLO','INTE',KLREO(I1234),MAX_STR_SPGP)
        CALL GETMEM('KLZ   ','ALLO','INTE',KLZ(I1234),LZ)
      END DO
*. Arrays for partitioning of Left vector = sigma
      NTTS = MXNTTS
      CALL GETMEM('LBT_L  ','ALLO','INTE',KLLBTL ,NTTS  )
      CALL GETMEM('LEBT_L ','ALLO','INTE',KLLEBTL,NTTS  )
      CALL GETMEM('I1BT_L ','ALLO','INTE',KLI1BTL,NTTS  )
      CALL GETMEM('IBT_L  ','ALLO','INTE',KLIBTL ,8*NTTS)
      CALL GETMEM('SCLF_L ','ALLO','REAL',KLSCLFCL,NTTS)
      CALL PART_CIV2(IDC,iWORK(KSBLTP),
     &               iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &               NOCTPA,NOCTPB,
     &               NSMST, LSCR1,
     &               iWORK(KSIOIO),ISMOST(1,ISSM),
     &               NBATCHL,
     &               iWORK(KLLBTL),iWORK(KLLEBTL),
     &               iWORK(KLI1BTL),iWORK(KLIBTL),
     &               0,ISIMSYM)
*. Number of BLOCKS
        NBLOCKL = IFRMR(IWORK(KLI1BTL),1,NBATCHL)
     &         + IFRMR(IWORK(KLLBTL),1,NBATCHL) - 1
*. Arrays for partitioning of Right  vector = C
      NTTS = MXNTTS
      CALL GETMEM('LBT_R  ','ALLO','INTE',KLLBTR ,NTTS  )
      CALL GETMEM('LEBT_R ','ALLO','INTE',KLLEBTR,NTTS  )
      CALL GETMEM('I1BT_R ','ALLO','INTE',KLI1BTR,NTTS  )
      CALL GETMEM('IBT_R  ','ALLO','INTE',KLIBTR ,8*NTTS)
      CALL GETMEM('SCLF_R ','ALLO','REAL',KLSCLFCR,NTTS)
      CALL PART_CIV2(IDC,iWORK(KCBLTP),
     &               iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &               NOCTPA,NOCTPB,
     &               NSMST, LSCR1,
     &               iWORK(KCIOIO),ISMOST(1,ICSM),
     &               NBATCHR,
     &               iWORK(KLLBTR),iWORK(KLLEBTR),
     &               iWORK(KLI1BTR),iWORK(KLIBTR),
     &               0,ISIMSYM)
*. Number of BLOCKS
        NBLOCKR = IFRMR(IWORK(KLI1BTR),1,NBATCHR)
     &         + IFRMR(IWORK(KLLBTR),1,NBATCHR) - 1
C?      WRITE(6,*) ' DENSI2T :NBLOCKR =',NBLOCKR


      IF(ICISTR.EQ.1) THEN
         WRITE(6,*) ' Sorry, ICISTR = 1 is out of fashion'
         WRITE(6,*) ' Switch to ICISTR = 2 - or reprogram '
*         STOP' DENSI2T : ICISTR = 1 in use '
         CALL SYSABENDMSG('lucia_util/densi2_lucia',
     &                    'Internal error',' ')
      ELSE IF(ICISTR.GE.2) THEN
        S2_TERM1 = 0.0D0
        CALL GASDN2_LUCIA(     I12,    RHO1,    RHO2,   RHO2S,   RHO2A,
     &                           L,       R,       L,       R,WORK(KC2),
     &                    iWORK(KCIOIO),iWORK(KSIOIO),
     &                    ISMOST(1,ICSM),ISMOST(1,ISSM),
     &                    iWORK(KCBLTP),iWORK(KSBLTP),NACOB,
     &                    iWORK(KNSTSO(IATP)),iWORK(KISTSO(IATP)),
     &                    iWORK(KNSTSO(IBTP)),iWORK(KISTSO(IBTP)),
     &                    NAEL,IATP,  NBEL,  IBTP,
     &                      IOCTPA,  IOCTPB,  NOCTPA,  NOCTPB,   NSMST,
     &                       NSMOB,   NSMSX,   NSMDX, MXPNGAS,  NOBPTS,
     &                      IOBPTS,    MAXK,    MAXI,   LSCR1,   LSCR1,
     &                    WORK(KCSCR),WORK(KSSCR),
     &                    SXSTSM,iWORK(KSTSTS),iWORK(KSTSTD),SXDXSX,
     &                    ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &                    iWORK(KI1),WORK(KXI1S),iWORK(KI2),WORK(KXI2S),
     &                    iWORK(KI3),WORK(KXI3S),iWORK(KI4),WORK(KXI4S),
     &                    WORK(KINSCR),MXPOBS,IPRDEN,WORK(KRHO1S),
     &                    LUL,LUR,PSSIGN,PSSIGN,
     &                    WORK(KRHO1P),WORK(KXNATO),
     &                    NBATCHL,
     &                    iWORK(KLLBTL),iWORK(KLLEBTL),
     &                    iWORK(KLI1BTL),iWORK(KLIBTL),
     &                    NBATCHR,
     &                    iWORK(KLLBTR),iWORK(KLLEBTR),
     &                    iWORK(KLI1BTR),iWORK(KLIBTR),
     &                    iWORK(KCONSPA),iWORK(KCONSPB),
     &                    WORK(KLSCLFCL),WORK(KLSCLFCR),
     &                    S2_TERM1, IUSE_PH,  IPHGAS,IDOSRHO1,   SRHO1,
     &                    IPACK)
*
        CALL GADSUM(RHO1,NACOB**2)
        IF(I12.EQ.2) THEN
          IF(IPACK) THEN
* If IPACK .EQ. .TRUE. then
C     Number of elements in symmetric and antisymmetric 2-body
C     density matrices are given in Nijkl.
            NIJ   = (NACOB*(NACOB+1))/2
            NIJKL = (NIJ*(NIJ+1))/2
            CALL GADSUM(RHO2S,NIJKL)
            CALL GADSUM(RHO2A,NIJKL)
          ELSE
            CALL GADSUM(RHO2,NACOB ** 2 *(NACOB**2+1)/2)
          END IF
        END IF
        IF(IDOSRHO1.EQ.1) THEN
          CALL GADSUM(SRHO1,NACOB ** 2)
        END IF
        CALL GADSUM_SCAL(S2_TERM1)
*
* CALL GASDN2_LUCIA --> 89
*
C     KLLBTR  KLLEBTR KLI1BTR KLIBTR
      END IF
C?    WRITE(6,*) ' Memcheck in densi2 after GASDN2'
C?    CALL MEMCHK
*
*
*. Add terms from hole-hole commutator
c      IF(IUSE_PH.EQ.1) THEN
c*. Overlap between left and right vector
c       XLR = INPRDD(L,R,LUR,LUL,1,-1)
c       CALL RHO1_HH(RHO1,XLR)
c      END IF

* Natural Orbitals
      CALL NATORB_LUCIA(RHO1,   NSMOB,  NTOOBS,  NACOBS,  NINOBS,IREOST,
     &                  WORK(KXNATO),WORK(KRHO1SM),WORK(KOCCSM),NACOB,
     &                  WORK(KRHO1P),IPRDEN)
*
      IF(IPRDEN.GE.5) THEN
        WRITE(6,*) ' One-electron density matrix '
        WRITE(6,*) ' ============================'
        CALL WRTMAT(RHO1,NTOOB,NTOOB,NTOOB,NTOOB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-electron density '
          CALL PRSYM(RHO2,NACOB**2)
        END IF
      END IF
*
      IF(I12.EQ.2) THEN
* <L!S**2|R>
        EXPS2 = S2_TERM1+0.25D0*DBLE(4*NAEL+(NAEL-NBEL)*(NAEL-NBEL-2))
        IF(IPRDEN.GT.0) THEN
          WRITE(6,*) ' Term 1 to S2 ', S2_TERM1
          WRITE(6,*) ' Expectation value of S2 ', EXPS2
        END IF
      ELSE
        EXPS2 = 0.0D0
      END IF
*
      IF(IDOSRHO1.EQ.1.AND.IPRDEN.GE.2) THEN
        WRITE(6,*) ' One-electron spindensity <0!E(aa) - E(bb)!0> '
        CALL WRTMAT(SRHO1,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF

*. Eliminate local memory
      CALL GETMEM('KSTSTS','FREE','INTE',KSTSTS,NSMST ** 2)
      CALL GETMEM('KSTSTD','FREE','INTE',KSTSTD,NSMST ** 2)
      CALL GETMEM('CONSPA','FREE','INTE',KCONSPA,NOCTPA**2)
      CALL GETMEM('CONSPB','FREE','INTE',KCONSPB,NOCTPB**2)
      CALL GETMEM('INSCR ','FREE','REAL',KINSCR,INTSCR)
      CALL GETMEM('SIOIO ','FREE','INTE',KSIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('CIOIO ','FREE','INTE',KCIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('I1    ','FREE','INTE',KI1,  LSCR3       )
      CALL GETMEM('I2    ','FREE','INTE',KI2,  LSCR3       )
      CALL GETMEM('I3    ','FREE','INTE',KI3,  LSCR3       )
      CALL GETMEM('I4    ','FREE','INTE',KI4,  LSCR3       )
      CALL GETMEM('XI1S  ','FREE','REAL',KXI1S,LSCR3       )
      CALL GETMEM('XI2S  ','FREE','REAL',KXI2S,LSCR3       )
      CALL GETMEM('XI3S  ','FREE','REAL',KXI3S,LSCR3       )
      CALL GETMEM('XI4S  ','FREE','REAL',KXI4S,LSCR3       )
      CALL GETMEM('SBLTP ','FREE','INTE',KSBLTP,NSMST)
      CALL GETMEM('CBLTP ','FREE','INTE',KCBLTP,NSMST)
      CALL GETMEM('RHO1S ','FREE','REAL',KRHO1S,NACOB ** 2)
      CALL GETMEM('RHO1P ','FREE','REAL',KRHO1P,NACOB*(NACOB+1)/2)
      CALL GETMEM('XNATO ','FREE','REAL',KXNATO,NACOB **2)
      CALL GETMEM('RHO1S ','FREE','REAL',KRHO1SM,NACOB ** 2)
      CALL GETMEM('RHO1S ','FREE','REAL',KXNATSM,NACOB ** 2)
      CALL GETMEM('RHO1S ','FREE','REAL',KOCCSM,NACOB )
      CALL GETMEM('KLZSCR','FREE','INTE',KLZSCR,LZSCR)
      DO K12 = 1, 1
        CALL GETMEM('KLOCS ','FREE','INTE',KLOCSTR(K12),MAX_STR_OC_BLK)
      END DO
      DO I1234 = 1, 2
        CALL GETMEM('KLREO ','FREE','INTE',KLREO(I1234),MAX_STR_SPGP)
        CALL GETMEM('KLZ   ','FREE','INTE',KLZ(I1234),LZ)
      END DO
      CALL GETMEM('LBT_L  ','FREE','INTE',KLLBTL ,NTTS  )
      CALL GETMEM('LEBT_L ','FREE','INTE',KLLEBTL,NTTS  )
      CALL GETMEM('I1BT_L ','FREE','INTE',KLI1BTL,NTTS  )
      CALL GETMEM('IBT_L  ','FREE','INTE',KLIBTL ,8*NTTS)
      CALL GETMEM('SCLF_L ','FREE','REAL',KLSCLFCL,NTTS)
      CALL GETMEM('LBT_R  ','FREE','INTE',KLLBTR ,NTTS  )
      CALL GETMEM('LEBT_R ','FREE','INTE',KLLEBTR,NTTS  )
      CALL GETMEM('I1BT_R ','FREE','INTE',KLI1BTR,NTTS  )
      CALL GETMEM('IBT_R  ','FREE','INTE',KLIBTR ,8*NTTS)
      CALL GETMEM('SCLF_R ','FREE','REAL',KLSCLFCR,NTTS)

      RETURN
      END
