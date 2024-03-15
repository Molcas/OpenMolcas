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
      use stdalloc, only: mma_allocate, mma_deallocate
      USE GLBBAS, only: VEC3
      use hidscr, only: ZSCR, ZOCSTR => OCSTR, REO, Z
      use strbas
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
#include "cands.fh"
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "crun.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "cprnt.fh"
#include "spinfo_lucia.fh"
*
      LOGICAL IPACK
#include "csmprd.fh"
#include "lucinp.fh"
#include "clunit.fh"
*. Scratch for string information
      INTEGER SXSTSM(1)
*. Specific input
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*),RHO2S(*),RHO2A(*),SRHO1(*)
      Integer, Allocatable:: CONSPA(:), CONSPB(:)
      Real*8, Allocatable:: INSCR(:)
      Integer, Allocatable:: STSTS(:), STSTD(:)
      Integer, Allocatable:: CIOIO(:), SIOIO(:)
      Integer, Allocatable:: CBLTP(:), SBLTP(:)
      Integer, Allocatable:: I1(:), I2(:), I3(:), I4(:)
      Real*8, Allocatable:: XI1S(:), XI2S(:), XI3S(:), XI4S(:)
      Integer, Allocatable:: LLBTL(:), LLBTR(:)
      Integer, Allocatable:: LLEBTL(:), LLEBTR(:)
      Integer, Allocatable:: LI1BTL(:), LI1BTR(:)
      Integer, Allocatable:: LIBTL(:), LIBTR(:)
      Real*8, Allocatable:: LSCLFCL(:), LSCLFCR(:)
      Integer, Allocatable:: SVST(:)
      Real*8, Allocatable:: RHO1S(:), RHO1P(:), XNATO(:), RHO1SM(:),
     &                       OCCSM(:)

*. Before I forget it :
*     IDUM = 0
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
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offsets for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)

* string sym, string sym => sx sym
* string sym, string sym => dx sym
      Call mma_allocate(STSTS,NSMST**2,Label='STSTS')
      Call mma_allocate(STSTD,NSMST**2,Label='STSTD')
      CALL STSTSM(STSTS,STSTD,NSMST)
*. connection matrices for supergroups
      Call mma_allocate(CONSPA,NOCTPA**2,Label='CONSPA')
      Call mma_allocate(CONSPB,NOCTPB**2,Label='CONSPB')
      CALL SPGRPCON(   IOCTPA,   NOCTPA,     NGAS,  MXPNGAS, NELFSPGP,
     &              CONSPA,IPRCIX)
      CALL SPGRPCON(   IOCTPB,   NOCTPB,     NGAS,  MXPNGAS, NELFSPGP,
     &              CONSPB,IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(NSTSO(IATP)%I,NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(NSTSO(IBTP)%I,NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(NSTSO(IATPM1)%I,NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(NSTSO(IATPM2)%I,NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(NSTSO(IBTPM1)%I,NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(NSTSO(IBTPM2)%I,NSMST*NOCTYP(IBTPM2),2)
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
* A 4 index block with four indices belonging OS class
      INTSCR = MXTSOB ** 4
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      Call mma_allocate(INSCR,INTSCR,Label='INSCR')
*
*. Arrays giving allowed type combinations '
      Call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
      Call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
*
      CALL IAIBCM(ISSPC,SIOIO)
      CALL IAIBCM(ISSPC,CIOIO)
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESCPH(CIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
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
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      IF (ENVIRO(1:6) .EQ. 'RASSCF') THEN
         KCSCR = LSCR12
      ELSE
         KCSCR = LSCR2
      END IF
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      Call mma_allocate(I1,LSCR3,Label='I1')
      Call mma_allocate(I2,LSCR3,Label='I2')
      Call mma_allocate(I3,LSCR3,Label='I3')
      Call mma_allocate(I4,LSCR3,Label='I4')
      Call mma_allocate(XI1S,LSCR3,Label='XI1S')
      Call mma_allocate(XI2S,LSCR3,Label='XI2S')
      Call mma_allocate(XI3S,LSCR3,Label='XI3S')
      Call mma_allocate(XI4S,LSCR3,Label='XI4S')
*. Arrays giving block type
      Call mma_allocate(SBLTP,NSMST,Label='SBLTP')
      Call mma_allocate(CBLTP,NSMST,Label='CBLTP')
*. Arrays for additional symmetry operation
c      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
c        Call mma_allocate(SVST,NSMST,Label='SVST')
c        CALL SIGVST(SVST,NSMST)
c      ELSE
         Call mma_allocate(SVST,1,Label='SVST')
c      END IF
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,SBLTP,SVST)
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,CBLTP,SVST)
      Call mma_deallocate(SVST)
* scratch space containing active one body
      CALL mma_allocate(RHO1S,NACOB ** 2,Label='RHO1S')
*. For natural orbitals
      CALL mma_allocate(RHO1P,NACOB*(NACOB+1)/2,Label='RHO1P')
      CALL mma_allocate(XNATO,NACOB **2,Label='XNATO')
*. Natural orbitals in symmetry blocks
      CALL mma_allocate(RHO1SM,NACOB ** 2,Label='RHO1SM')
      CALL mma_allocate(OCCSM,NACOB,Label='OCCSM')
*
*. Space for one block of string occupations and two arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      call mma_allocate(ZSCR,lZSCR,Label='ZSCR')
      K12=1
      call mma_allocate(ZOCSTR,MAX_STR_OC_BLK,K12,Label='ZOCSTR')
      I1234=2
      Call mma_allocate(REO,MAX_STR_SPGP,I1234,Label='REO')
      CALL mma_allocate(Z,LZ,I1234,Label='Z')
*. Arrays for partitioning of Left vector = sigma
      NTTS = MXNTTS
      Call mma_allocate(LLBTL,NTTS,Label='LLBTL')
      Call mma_allocate(LLEBTL,NTTS,Label='LLEBTL')
      Call mma_allocate(LI1BTL,NTTS,Label='LI1BTL')
      Call mma_allocate(LIBTL,8*NTTS,Label='LIBTL')
      Call mma_allocate(LSCLFCL,NTTS,Label='LSCLFCL')
      CALL PART_CIV2(IDC,SBLTP,
     &               NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &               NOCTPA,NOCTPB,
     &               NSMST, LSCR1,
     &               SIOIO,ISMOST(1,ISSM),
     &               NBATCHL,
     &               LLBTL,LLEBTL,
     &               LI1BTL,LIBTL,
     &               0,ISIMSYM)
*. Arrays for partitioning of Right  vector = C
      NTTS = MXNTTS
      Call mma_allocate(LLBTR,NTTS,Label='LLBTR')
      Call mma_allocate(LLEBTR,NTTS,Label='LLEBTR')
      Call mma_allocate(LI1BTR,NTTS,Label='LI1BTR')
      Call mma_allocate(LIBTR,8*NTTS,Label='LIBTR')
      Call mma_allocate(LSCLFCR,NTTS,Label='LSCLFCR')
      CALL PART_CIV2(IDC,CBLTP,
     &               NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &               NOCTPA,NOCTPB,
     &               NSMST, LSCR1,
     &               CIOIO,ISMOST(1,ICSM),
     &               NBATCHR,
     &               LLBTR,LLEBTR,
     &               LI1BTR,LIBTR,
     &               0,ISIMSYM)

      IF(ICISTR.EQ.1) THEN
         WRITE(6,*) ' Sorry, ICISTR = 1 is out of fashion'
         WRITE(6,*) ' Switch to ICISTR = 2 - or reprogram '
*         STOP' DENSI2T : ICISTR = 1 in use '
         CALL SYSABENDMSG('lucia_util/densi2_lucia',
     &                    'Internal error',' ')
      ELSE IF(ICISTR.GE.2) THEN
        S2_TERM1 = 0.0D0
        CALL GASDN2_LUCIA(     I12,    RHO1,    RHO2,   RHO2S,   RHO2A,
     &                           L,       R,       L,     R,VEC3,
     &                    CIOIO,SIOIO,
     &                    ISMOST(1,ICSM),ISMOST(1,ISSM),
     &                    CBLTP,SBLTP,NACOB,
     &                    NSTSO(IATP)%I,ISTSO(IATP)%I,
     &                    NSTSO(IBTP)%I,ISTSO(IBTP)%I,
     &                    NAEL,IATP,  NBEL,  IBTP,
     &                      IOCTPA,  IOCTPB,  NOCTPA,  NOCTPB,   NSMST,
     &                       NSMOB,   NSMSX,   NSMDX, MXPNGAS,  NOBPTS,
     &                      IOBPTS,    MAXK,    MAXI,   LSCR1,   LSCR1,
     &                    VEC3(1+KCSCR),VEC3,
     &                    SXSTSM,STSTS,STSTD,SXDXSX,
     &                    ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &                    I1,XI1S,I2,XI2S,
     &                    I3,XI3S,I4,XI4S,
     &                    INSCR,MXPOBS,IPRDEN,RHO1S,
     &                    LUL,LUR,PSSIGN,PSSIGN,
     &                    RHO1P,XNATO,
     &                    NBATCHL,
     &                    LLBTL,LLEBTL,
     &                    LI1BTL,LIBTL,
     &                    NBATCHR,
     &                    LLBTR,LLEBTR,
     &                    LI1BTR,LIBTR,
     &                    CONSPA,CONSPB,
     &                    LSCLFCL,LSCLFCR,
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
C     LBTR  LLEBTR LI1BTR LIBTR
      END IF
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
     &                  XNATO,RHO1SM,OCCSM,NACOB,
     &                  RHO1P,IPRDEN)
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
      Call mma_deallocate(STSTS)
      Call mma_deallocate(STSTD)
      Call mma_deallocate(CONSPA)
      Call mma_deallocate(CONSPB)
      Call mma_deallocate(INSCR)
      Call mma_deallocate(SIOIO)
      Call mma_deallocate(CIOIO)
      Call mma_deallocate(I1)
      Call mma_deallocate(I2)
      Call mma_deallocate(I3)
      Call mma_deallocate(I4)
      Call mma_deallocate(XI1S)
      Call mma_deallocate(XI2S)
      Call mma_deallocate(XI3S)
      Call mma_deallocate(XI4S)
      Call mma_deallocate(SBLTP)
      Call mma_deallocate(CBLTP)
      Call mma_deallocate(RHO1S)
      Call mma_deallocate(RHO1P)
      Call mma_deallocate(XNATO)
      Call mma_deallocate(RHO1SM)
      Call mma_deallocate(OCCSM)
      Call mma_deallocate(ZSCR)
      Call mma_deallocate(ZOCSTR)
      Call mma_deallocate(REO)
      Call mma_deallocate(Z)
      Call mma_deallocate(LLBTL)
      Call mma_deallocate(LLEBTL)
      Call mma_deallocate(LI1BTL)
      Call mma_deallocate(LIBTL)
      Call mma_deallocate(LSCLFCL)
      Call mma_deallocate(LLBTR)
      Call mma_deallocate(LLEBTR)
      Call mma_deallocate(LI1BTR)
      Call mma_deallocate(LIBTR)
      Call mma_deallocate(LSCLFCR)

      RETURN
      END
