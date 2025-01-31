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

subroutine DENSI2_LUCIA(I12,RHO1,RHO2,RHO2S,RHO2A,L,R,LUL,LUR,EXPS2,IDOSRHO1,SRHO1,IPACK)
! Density matrices between L and R
!
! I12 = 1 => only one-body density
! I12 = 2 => one- and two-body density matrices
!
! Jeppe Olsen,      Oct 94
! GAS modifications Aug 95
! Two body density added, '96
!
! Table-Block driven, June 97
! Spin density added, Jan. 99
!
! Jesper Wisborg Krogh
! Allowing to symmetry pack on the fly, Sept. 2003
!
! Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
! ijkl = ij*(ij-1)/2+kl, ij >= kl
!
! Two-body symmetric density stored in rho2s
! Two-body anti-symmetric density stored in rho2a
!
! If the two-body density matrix is calculated, then also the
! expectation value of the spin is evaluated.
! The latter is realized as
! S**2
!      = S+S- + Sz(Sz-1)
!      = -Sum(ij) a+i alpha a+j beta a i beta a j alpha + Nalpha +
!        1/2(N alpha - N beta))(1/2(N alpha - Nbeta) - 1)
! If IDOSRHO1 = 1, spin density is also calculated
!
! =====
! Input
! =====
!
!.Definition of L and R is picked up from CANDS
! with L being S and  R being C
!

use stdalloc, only: mma_allocate, mma_deallocate
use GLBBAS, only: VEC3
use hidscr, only: ZSCR, ZOCSTR => OCSTR, REO, Z
use strbas, only: NSTSO, ISTSO
use CandS, only: ICSM, ISSM, ISSPC
use Constants, only: Zero
use lucia_data, only: NGAS, IPHGAS
use lucia_data, only: MXSB, MXSOOB, MXNTTS, ISMOST, XISPSM
use lucia_data, only: IPRDEN, IPRCIX
use lucia_data, only: ENVIRO, ICISTR, IADVICE, ISIMSYM, IUSE_PH, LCSBLK, MXINKA
use lucia_data, only: IREFSM, PSSIGN, IDC
use lucia_data, only: MXNSTR, IBSPGPFTP, MAX_STR_OC_BLK, MAX_STR_SPGP, MNHL, NELFSPGP, NHLFSPGP, NSTFSMSPGP
use lucia_data, only: NSMOB
use lucia_data, only: NACOB, MXTSOB, NOCOB, IOBPTS, IREOST, NACOBS, NINOBS, NOBPTS, NTOOB, NTOOBS
use lucia_data, only: NOCTYP
use lucia_data, only: NELEC
use lucia_data, only: MXPOBS, MXPNGAS, MXPNSMST
use csm_data, only: NSMST, NSMDX, NSMSX
use csm_data, only: ADSXA, ASXAD, SXDXSX

implicit none
! Specific input
integer I12, LUL, LUR, IDOSRHO1
logical IPACK
real*8 L(*), R(*)
! Output
real*8 RHO1(*), RHO2(*), RHO2S(*), RHO2A(*), SRHO1(*)
real*8 EXPS2
integer, allocatable :: CONSPA(:), CONSPB(:)
real*8, allocatable :: INSCR(:)
integer, allocatable :: STSTS(:), STSTD(:)
integer, allocatable :: CIOIO(:), SIOIO(:)
integer, allocatable :: CBLTP(:), SBLTP(:)
integer, allocatable :: I1(:), I2(:), I3(:), I4(:)
real*8, allocatable :: XI1S(:), XI2S(:), XI3S(:), XI4S(:)
integer, allocatable :: LLBTL(:), LLBTR(:)
integer, allocatable :: LLEBTL(:), LLEBTR(:)
integer, allocatable :: LI1BTL(:), LI1BTR(:)
integer, allocatable :: LIBTL(:), LIBTR(:)
real*8, allocatable :: LSCLFCL(:), LSCLFCR(:)
integer, allocatable :: SVST(:)
real*8, allocatable :: RHO1S(:), RHO1P(:), XNATO(:), RHO1SM(:), OCCSM(:)
! Scratch for string information
integer SXSTSM(1)
integer NIJ, NIJKL, IATP, IBTP, IATPM1, IBTPM1, IATPM2, IBTPM2, NOCTPA, NOCTPB, IOCTPA, IOCTPB, NAEL, NBEL, MAXA0, MAXB0, MXSTBL0, &
        MAXA, MAXA1, MAXB, MAXB1, MXSTBL, MAXI, MAXK, IOBTP, IOBSM, LSCR1, INTSCR, MXCJ, MXCIJA, MXCIJB, MXCIJAB, MXSXBL, LSCR2, &
        LSCR12, KCSCR, MAXIK, LSCR3, LZSCR, LZ, K12, I1234, NTTS
integer MXADKBLK, MXADKBLK_AS, MXCJ_ALLSYM, MX_NSPII, NBATCHL, NBATCHR
integer, external :: IMNMX
real*8 S2_TERM1

! Before I forget it :
!IDUM = 0
!call MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
call SETVEC(RHO1,ZERO,NACOB**2)
if (I12 == 2) then
  if (IPACK) then
    ! If IPACK == .TRUE. then
    ! Number of elements in symmetric and antisymmetric 2-body
    ! density matrices are given in Nijkl.
    NIJ = (NACOB*(NACOB+1))/2
    NIJKL = (NIJ*(NIJ+1))/2
    call SETVEC(RHO2S,ZERO,NIJKL)
    call SETVEC(RHO2A,ZERO,NIJKL)
  else
    call SETVEC(RHO2,ZERO,NACOB**2*(NACOB**2+1)/2)
  end if
end if

if (IDOSRHO1 == 1) call SETVEC(SRHO1,ZERO,NACOB**2)

! Info for this internal space
! type of alpha and beta strings
IATP = 1
IBTP = 2
! alpha and beta strings with an electron removed
IATPM1 = 3
IBTPM1 = 4
! alpha and beta strings with two electrons removed
IATPM2 = 5
IBTPM2 = 6
! Number of supergroups
NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Offsets for supergroups
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)

! string sym, string sym => sx sym
! string sym, string sym => dx sym
call mma_allocate(STSTS,NSMST**2,Label='STSTS')
call mma_allocate(STSTD,NSMST**2,Label='STSTD')
call STSTSM(STSTS,STSTD,NSMST)
! connection matrices for supergroups
call mma_allocate(CONSPA,NOCTPA**2,Label='CONSPA')
call mma_allocate(CONSPB,NOCTPB**2,Label='CONSPB')
call SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,CONSPA,IPRCIX)
call SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,CONSPB,IPRCIX)
! Largest block of strings in zero order space
MAXA0 = IMNMX(NSTSO(IATP)%I,NSMST*NOCTYP(IATP),2)
MAXB0 = IMNMX(NSTSO(IBTP)%I,NSMST*NOCTYP(IBTP),2)
MXSTBL0 = MXNSTR
! Largest number of strings of given symmetry and type
MAXA = 0
if (NAEL >= 1) then
  MAXA1 = IMNMX(NSTSO(IATPM1)%I,NSMST*NOCTYP(IATPM1),2)
  MAXA = max(MAXA,MAXA1)
  if (NAEL >= 2) then
    MAXA1 = IMNMX(NSTSO(IATPM2)%I,NSMST*NOCTYP(IATPM2),2)
    MAXA = max(MAXA,MAXA1)
  end if
end if
MAXB = 0
if (NBEL >= 1) then
  MAXB1 = IMNMX(NSTSO(IBTPM1)%I,NSMST*NOCTYP(IBTPM1),2)
  MAXB = max(MAXB,MAXB1)
  if (NBEL >= 2) then
    MAXB1 = IMNMX(NSTSO(IBTPM2)%I,NSMST*NOCTYP(IBTPM2),2)
    MAXB = max(MAXB,MAXB1)
  end if
end if
MAXA = max(MAXA,MAXA0)
MAXB = max(MAXB,MAXB0)
MXSTBL = max(MAXA,MAXB)
if (IPRDEN >= 2) write(6,*) ' Largest block of strings with given symmetry and type',MXSTBL
! Largest number of resolution strings and spectator strings
! that can be treated simultaneously
! replace with MXINKA !!!
MAXI = min(MXINKA,MXSTBL)
MAXK = min(MXINKA,MXSTBL)
!write(6,*) ' DENSI2 : MAXI MAXK ',MAXI,MAXK
! Largest active orbital block belonging to given type and symmetry
MXTSOB = 0
do IOBTP=1,NGAS
  do IOBSM=1,NSMOB
    MXTSOB = max(MXTSOB,NOBPTS(IOBTP,IOBSM))
  end do
end do
! Local scratch arrays for blocks of C and sigma
if (IPRDEN >= 2) write(6,*) ' DENSI2 : MXSB MXTSOB MXSOOB ',MXSB,MXTSOB,MXSOOB
!if (ISIMSYM /= 1) THEN
LSCR1 = MXSOOB
!else
!  LSCR1 = MXSOOB_AS
!end if
LSCR1 = max(LSCR1,LCSBLK)
! JESPER: Should reduce I/O
if (ENVIRO(1:6) == 'RASSCF') then
  LSCR1 = max(int(XISPSM(IREFSM,1)),MXSOOB)
  if (PSSIGN /= 0.0d0) LSCR1 = int(2.0d0*XISPSM(IREFSM,1))
end if
if (IPRDEN >= 2) write(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
! SCRATCH space for block of two-electron density matrix
! A 4 index block with four indices belonging OS class
INTSCR = MXTSOB**4
if (IPRDEN >= 2) write(6,*) ' Density scratch space ',INTSCR
call mma_allocate(INSCR,INTSCR,Label='INSCR')

! Arrays giving allowed type combinations
call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')

call IAIBCM(ISSPC,SIOIO)
call IAIBCM(ISSPC,CIOIO)
! Scratch space for CJKAIB resolution matrices
call MXRESCPH(CIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSMST,NSTFSMSPGP,MXPNSMST,NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,NELFSPGP,MXCJ, &
              MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
if (IPRDEN >= 2) write(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
LSCR2 = max(MXCJ,MXCIJA,MXCIJB)
if (IPRDEN >= 2) write(6,*) ' Space for resolution matrices ',LSCR2
LSCR12 = max(LSCR1,2*LSCR2)
if (ENVIRO(1:6) == 'RASSCF') LSCR12 = max(LSCR1,LSCR2)
! It is assumed that the third block already has been allocated, so
if (IPRCIX >= 2) write(6,*) ' Space for resolution matrices ',LSCR12
if (ENVIRO(1:6) == 'RASSCF') then
  KCSCR = LSCR12
else
  KCSCR = LSCR2
end if

! Space for annihilation/creation mappings
MAXIK = max(MAXI,MAXK)
LSCR3 = max(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
call mma_allocate(I1,LSCR3,Label='I1')
call mma_allocate(I2,LSCR3,Label='I2')
call mma_allocate(I3,LSCR3,Label='I3')
call mma_allocate(I4,LSCR3,Label='I4')
call mma_allocate(XI1S,LSCR3,Label='XI1S')
call mma_allocate(XI2S,LSCR3,Label='XI2S')
call mma_allocate(XI3S,LSCR3,Label='XI3S')
call mma_allocate(XI4S,LSCR3,Label='XI4S')
! Arrays giving block type
call mma_allocate(SBLTP,NSMST,Label='SBLTP')
call mma_allocate(CBLTP,NSMST,Label='CBLTP')
! Arrays for additional symmetry operation
!if ((IDC == 3) .or. (IDC == 4)) then
!  call mma_allocate(SVST,NSMST,Label='SVST')
!  call SIGVST(SVST,NSMST)
!else
call mma_allocate(SVST,1,Label='SVST')
!end if
call ZBLTP(ISMOST(1,ISSM),NSMST,IDC,SBLTP,SVST)
call ZBLTP(ISMOST(1,ICSM),NSMST,IDC,CBLTP,SVST)
call mma_deallocate(SVST)
! scratch space containing active one body
call mma_allocate(RHO1S,NACOB**2,Label='RHO1S')
! For natural orbitals
call mma_allocate(RHO1P,NACOB*(NACOB+1)/2,Label='RHO1P')
call mma_allocate(XNATO,NACOB**2,Label='XNATO')
! Natural orbitals in symmetry blocks
call mma_allocate(RHO1SM,NACOB**2,Label='RHO1SM')
call mma_allocate(OCCSM,NACOB,Label='OCCSM')

! Space for one block of string occupations and two arrays of
! reordering arrays
LZSCR = (max(NAEL,NBEL)+3)*(NOCOB+1)+2*NOCOB
LZ = (max(NAEL,NBEL)+2)*NOCOB
call mma_allocate(ZSCR,lZSCR,Label='ZSCR')
K12 = 1
call mma_allocate(ZOCSTR,MAX_STR_OC_BLK,K12,Label='ZOCSTR')
I1234 = 2
call mma_allocate(REO,MAX_STR_SPGP,I1234,Label='REO')
call mma_allocate(Z,LZ,I1234,Label='Z')
! Arrays for partitioning of Left vector = sigma
NTTS = MXNTTS
call mma_allocate(LLBTL,NTTS,Label='LLBTL')
call mma_allocate(LLEBTL,NTTS,Label='LLEBTL')
call mma_allocate(LI1BTL,NTTS,Label='LI1BTL')
call mma_allocate(LIBTL,8*NTTS,Label='LIBTL')
call mma_allocate(LSCLFCL,NTTS,Label='LSCLFCL')
call PART_CIV2(IDC,SBLTP,NSTSO(IATP)%I,NSTSO(IBTP)%I,NOCTPA,NOCTPB,NSMST,LSCR1,SIOIO,ISMOST(1,ISSM),NBATCHL,LLBTL,LLEBTL,LI1BTL, &
               LIBTL,0,ISIMSYM)
! Arrays for partitioning of Right  vector = C
NTTS = MXNTTS
call mma_allocate(LLBTR,NTTS,Label='LLBTR')
call mma_allocate(LLEBTR,NTTS,Label='LLEBTR')
call mma_allocate(LI1BTR,NTTS,Label='LI1BTR')
call mma_allocate(LIBTR,8*NTTS,Label='LIBTR')
call mma_allocate(LSCLFCR,NTTS,Label='LSCLFCR')
call PART_CIV2(IDC,CBLTP,NSTSO(IATP)%I,NSTSO(IBTP)%I,NOCTPA,NOCTPB,NSMST,LSCR1,CIOIO,ISMOST(1,ICSM),NBATCHR,LLBTR,LLEBTR,LI1BTR, &
               LIBTR,0,ISIMSYM)

if (ICISTR == 1) then
  write(6,*) ' Sorry, ICISTR = 1 is out of fashion'
  write(6,*) ' Switch to ICISTR = 2 - or reprogram'
  !stop ' DENSI2T : ICISTR = 1 in use'
  call SYSABENDMSG('lucia_util/densi2_lucia','Internal error','')
else if (ICISTR >= 2) then
  S2_TERM1 = 0.0d0
  call GASDN2_LUCIA(I12,RHO1,RHO2,RHO2S,RHO2A,L,R,L,R,VEC3,CIOIO,SIOIO,ISMOST(1,ICSM),ISMOST(1,ISSM),CBLTP,SBLTP,NACOB, &
                    NSTSO(IATP)%I,ISTSO(IATP)%I,NSTSO(IBTP)%I,ISTSO(IBTP)%I,NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSMST, &
                    NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,MAXK,MAXI,LSCR1,LSCR1,VEC3(1+KCSCR),VEC3,SXSTSM,STSTS,STSTD,SXDXSX, &
                    ADSXA,ASXAD,NGAS,NELFSPGP,IDC,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,INSCR,MXPOBS,IPRDEN,RHO1S,LUL,LUR,PSSIGN,PSSIGN, &
                    RHO1P,XNATO,NBATCHL,LLBTL,LLEBTL,LI1BTL,LIBTL,NBATCHR,LLBTR,LLEBTR,LI1BTR,LIBTR,CONSPA,CONSPB,LSCLFCL,LSCLFCR, &
                    S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IPACK)

  call GADSUM(RHO1,NACOB**2)
  if (I12 == 2) then
    if (IPACK) then
      ! If IPACK == .TRUE. then
      ! Number of elements in symmetric and antisymmetric 2-body
      ! density matrices are given in Nijkl.
      NIJ = (NACOB*(NACOB+1))/2
      NIJKL = (NIJ*(NIJ+1))/2
      call GADSUM(RHO2S,NIJKL)
      call GADSUM(RHO2A,NIJKL)
    else
      call GADSUM(RHO2,NACOB**2*(NACOB**2+1)/2)
    end if
  end if
  if (IDOSRHO1 == 1) call GADSUM(SRHO1,NACOB**2)
  call GADSUM_SCAL(S2_TERM1)

  ! CALL GASDN2_LUCIA --> 89
  !
  ! LBTR  LLEBTR LI1BTR LIBTR
end if

! Add terms from hole-hole commutator
!if (IUSE_PH == 1) then
!  !Overlap between left and right vector
!  XLR = INPRDD(L,R,LUR,LUL,1,-1)
!  call RHO1_HH(RHO1,XLR)
!end if

! Natural Orbitals
call NATORB_LUCIA(RHO1,NSMOB,NTOOBS,NACOBS,NINOBS,IREOST,XNATO,RHO1SM,OCCSM,NACOB,RHO1P,IPRDEN)

if (IPRDEN >= 5) then
  write(6,*) ' One-electron density matrix'
  write(6,*) ' ==========================='
  call WRTMAT(RHO1,NTOOB,NTOOB,NTOOB,NTOOB)
  if (I12 == 2) then
    write(6,*) ' Two-electron density'
    call PRSYM(RHO2,NACOB**2)
  end if
end if

if (I12 == 2) then
  ! <L!S**2|R>
  EXPS2 = S2_TERM1+0.25d0*dble(4*NAEL+(NAEL-NBEL)*(NAEL-NBEL-2))
  if (IPRDEN > 0) then
    write(6,*) ' Term 1 to S2 ',S2_TERM1
    write(6,*) ' Expectation value of S2 ',EXPS2
  end if
else
  EXPS2 = 0.0d0
end if

if ((IDOSRHO1 == 1) .and. (IPRDEN >= 2)) then
  write(6,*) ' One-electron spindensity <0!E(aa) - E(bb)!0>'
  call WRTMAT(SRHO1,NTOOB,NTOOB,NTOOB,NTOOB)
end if

! Eliminate local memory
call mma_deallocate(STSTS)
call mma_deallocate(STSTD)
call mma_deallocate(CONSPA)
call mma_deallocate(CONSPB)
call mma_deallocate(INSCR)
call mma_deallocate(SIOIO)
call mma_deallocate(CIOIO)
call mma_deallocate(I1)
call mma_deallocate(I2)
call mma_deallocate(I3)
call mma_deallocate(I4)
call mma_deallocate(XI1S)
call mma_deallocate(XI2S)
call mma_deallocate(XI3S)
call mma_deallocate(XI4S)
call mma_deallocate(SBLTP)
call mma_deallocate(CBLTP)
call mma_deallocate(RHO1S)
call mma_deallocate(RHO1P)
call mma_deallocate(XNATO)
call mma_deallocate(RHO1SM)
call mma_deallocate(OCCSM)
call mma_deallocate(ZSCR)
call mma_deallocate(ZOCSTR)
call mma_deallocate(REO)
call mma_deallocate(Z)
call mma_deallocate(LLBTL)
call mma_deallocate(LLEBTL)
call mma_deallocate(LI1BTL)
call mma_deallocate(LIBTL)
call mma_deallocate(LSCLFCL)
call mma_deallocate(LLBTR)
call mma_deallocate(LLEBTR)
call mma_deallocate(LI1BTR)
call mma_deallocate(LIBTR)
call mma_deallocate(LSCLFCR)

end subroutine DENSI2_LUCIA
