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

subroutine DENSI2(I12,RHO1,RHO2,RHO2S,RHO2A,L,R,LUL,LUR,EXPS2,IDOSRHO1,SRHO1,IPACK)
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
! ijkl = nTri_Elem(ij)+kl, ij >= kl
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
! Definition of L and R is picked up from CANDS
! with L being S and R being C
!

use Index_Functions, only: nTri_Elem
use CandS, only: ICSM, ISSM, ISSPC
use lucia_data, only: ENVIRO, IADVICE, IBSPGPFTP, ICISTR, IDC, IOBPTS, IPHGAS, IPRCIX, IPRDEN, IREFSM, IREOST, LCSBLK, &
                      MAX_STR_OC_BLK, MAX_STR_SPGP, MNHL, MXINKA, MXNSTR, MXNTTS, MXPNGAS, MXPNSMST, MXSOOB, MXTSOB, NACOB, &
                      NACOBS, NELEC, NELFSPGP, NGAS, NHLFSPGP, NINOBS, NIRREP, NOBPTS, NOCOB, NOCTYP, NSMOB, NSTFSMSPGP, NSTSO, &
                      NTOOB, NTOOBS, OCSTR, PSSIGN, REO, VEC3, XISPSM, Z, ZSCR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: I12, LUL, LUR, IDOSRHO1
real(kind=wp), intent(inout) :: RHO1(*), RHO2(*), RHO2S(*), RHO2A(*), L(*), R(*), SRHO1(*)
real(kind=wp), intent(out) :: EXPS2
logical(kind=iwp), intent(in) :: IPACK
integer(kind=iwp) :: I1234, IATP, IATPM1, IATPM2, IBTP, IBTPM1, IBTPM2, INTSCR, IOCTPA, IOCTPB, K12, KCSCR, LSCR1, LSCR12, LSCR2, &
                     LSCR3, LZ, LZSCR, MAXA, MAXA0, MAXA1, MAXB, MAXB0, MAXB1, MAXI, MAXIK, MAXK, MX_NSPII, MXADKBLK, MXADKBLK_AS, &
                     MXCIJA, MXCIJAB, MXCIJB, MXCJ, MXCJ_ALLSYM, MXSTBL, MXSTBL0, MXSXBL, NAEL, NBATCHL, NBATCHR, NBEL, NIJ, &
                     NIJKL, NOCTPA, NOCTPB, NTTS
real(kind=wp) :: S2_TERM1
integer(kind=iwp), allocatable :: CBLTP(:), CIOIO(:), CONSPA(:), CONSPB(:), I1(:), I2(:), I3(:), I4(:), LI1BTL(:), LI1BTR(:), &
                                  LIBTL(:), LIBTR(:), LLBTL(:), LLBTR(:), LLEBTL(:), LLEBTR(:), SBLTP(:), SIOIO(:), SVST(:)
real(kind=wp), allocatable :: INSCR(:), LSCLFCL(:), LSCLFCR(:), OCCSM(:), RHO1P(:), RHO1S(:), RHO1SM(:), XI1S(:), XI2S(:), &
                              XI3S(:), XI4S(:), XNATO(:)
integer(kind=iwp), external :: IMNMX

! Before I forget it :
!IDUM = 0
!call MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
RHO1(1:NACOB**2) = Zero
if (I12 == 2) then
  if (IPACK) then
    ! If IPACK == .TRUE. then
    ! Number of elements in symmetric and antisymmetric 2-body
    ! density matrices are given in Nijkl.
    NIJ = nTri_Elem(NACOB)
    NIJKL = nTri_Elem(NIJ)
    RHO2S(1:NIJKL) = Zero
    RHO2A(1:NIJKL) = Zero
  else
    RHO2(1:nTri_Elem(NACOB**2)) = Zero
  end if
end if

if (IDOSRHO1 == 1) SRHO1(1:NACOB**2) = Zero

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

! connection matrices for supergroups
call mma_allocate(CONSPA,NOCTPA**2,Label='CONSPA')
call mma_allocate(CONSPB,NOCTPB**2,Label='CONSPB')
call SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,CONSPA)
call SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,CONSPB)
! Largest block of strings in zero order space
MAXA0 = IMNMX(NSTSO(IATP)%A,NIRREP*NOCTYP(IATP),2)
MAXB0 = IMNMX(NSTSO(IBTP)%A,NIRREP*NOCTYP(IBTP),2)
MXSTBL0 = MXNSTR
! Largest number of strings of given symmetry and type
MAXA = 0
if (NAEL >= 1) then
  MAXA1 = IMNMX(NSTSO(IATPM1)%A,NIRREP*NOCTYP(IATPM1),2)
  MAXA = max(MAXA,MAXA1)
  if (NAEL >= 2) then
    MAXA1 = IMNMX(NSTSO(IATPM2)%A,NIRREP*NOCTYP(IATPM2),2)
    MAXA = max(MAXA,MAXA1)
  end if
end if
MAXB = 0
if (NBEL >= 1) then
  MAXB1 = IMNMX(NSTSO(IBTPM1)%A,NIRREP*NOCTYP(IBTPM1),2)
  MAXB = max(MAXB,MAXB1)
  if (NBEL >= 2) then
    MAXB1 = IMNMX(NSTSO(IBTPM2)%A,NIRREP*NOCTYP(IBTPM2),2)
    MAXB = max(MAXB,MAXB1)
  end if
end if
MAXA = max(MAXA,MAXA0)
MAXB = max(MAXB,MAXB0)
MXSTBL = max(MAXA,MAXB)
if (IPRDEN >= 2) write(u6,*) ' Largest block of strings with given symmetry and type',MXSTBL
! Largest number of resolution strings and spectator strings
! that can be treated simultaneously
! replace with MXINKA !!!
MAXI = min(MXINKA,MXSTBL)
MAXK = min(MXINKA,MXSTBL)
!write(u6,*) ' DENSI2 : MAXI MAXK ',MAXI,MAXK
! Largest active orbital block belonging to given type and symmetry
MXTSOB = max(0,maxval(NOBPTS(1:NGAS,1:NSMOB)))
! Local scratch arrays for blocks of C and sigma
if (IPRDEN >= 2) write(u6,*) ' DENSI2 : MXTSOB MXSOOB ',MXTSOB,MXSOOB
!if (ISIMSYM /= 1) THEN
LSCR1 = MXSOOB
!else
!  LSCR1 = MXSOOB_AS
!end if
LSCR1 = max(LSCR1,LCSBLK)
! JESPER: Should reduce I/O
if (ENVIRO == 'RASSCF') then
  LSCR1 = max(int(XISPSM(IREFSM,1)),MXSOOB)
  if (PSSIGN /= Zero) LSCR1 = int(Two*XISPSM(IREFSM,1))
end if
if (IPRDEN >= 2) write(u6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
! SCRATCH space for block of two-electron density matrix
! A 4 index block with four indices belonging OS class
INTSCR = MXTSOB**4
if (IPRDEN >= 2) write(u6,*) ' Density scratch space ',INTSCR
call mma_allocate(INSCR,INTSCR,Label='INSCR')

! Arrays giving allowed type combinations
call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')

call IAIBCM(ISSPC,SIOIO)
call IAIBCM(ISSPC,CIOIO)
! Scratch space for CJKAIB resolution matrices
call MXRESCPH(CIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NIRREP,NSTFSMSPGP,MXPNSMST,NSMOB,MXPNGAS,NGAS,NOBPTS,MAXK,NELFSPGP,MXCJ,MXCIJA, &
              MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
if (IPRDEN >= 2) write(u6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
LSCR2 = max(MXCJ,MXCIJA,MXCIJB)
if (IPRDEN >= 2) write(u6,*) ' Space for resolution matrices ',LSCR2
LSCR12 = max(LSCR1,2*LSCR2)
if (ENVIRO == 'RASSCF') LSCR12 = max(LSCR1,LSCR2)
! It is assumed that the third block already has been allocated, so
if (IPRCIX >= 2) write(u6,*) ' Space for resolution matrices ',LSCR12
if (ENVIRO == 'RASSCF') then
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
call mma_allocate(SBLTP,NIRREP,Label='SBLTP')
call mma_allocate(CBLTP,NIRREP,Label='CBLTP')
! Arrays for additional symmetry operation
!if ((IDC == 3) .or. (IDC == 4)) then
!  call mma_allocate(SVST,NIRREP,Label='SVST')
!  call SIGVST(SVST,NIRREP)
!else
call mma_allocate(SVST,1,Label='SVST')
!end if
call ZBLTP(ISSM,NIRREP,IDC,SBLTP,SVST)
call ZBLTP(ICSM,NIRREP,IDC,CBLTP,SVST)
call mma_deallocate(SVST)
! scratch space containing active one body
call mma_allocate(RHO1S,NACOB**2,Label='RHO1S')
! For natural orbitals
call mma_allocate(RHO1P,nTri_Elem(NACOB),Label='RHO1P')
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
call mma_allocate(OCSTR,MAX_STR_OC_BLK,K12,Label='OCSTR')
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
call PART_CIV2(IDC,NSTSO(IATP)%A,NSTSO(IBTP)%A,NOCTPA,NOCTPB,NIRREP,SIOIO,ISSM,NBATCHL,LLBTL,LLEBTL,LI1BTL,LIBTL,0)
! Arrays for partitioning of Right  vector = C
NTTS = MXNTTS
call mma_allocate(LLBTR,NTTS,Label='LLBTR')
call mma_allocate(LLEBTR,NTTS,Label='LLEBTR')
call mma_allocate(LI1BTR,NTTS,Label='LI1BTR')
call mma_allocate(LIBTR,8*NTTS,Label='LIBTR')
call mma_allocate(LSCLFCR,NTTS,Label='LSCLFCR')
call PART_CIV2(IDC,NSTSO(IATP)%A,NSTSO(IBTP)%A,NOCTPA,NOCTPB,NIRREP,CIOIO,ICSM,NBATCHR,LLBTR,LLEBTR,LI1BTR,LIBTR,0)

if (ICISTR == 1) then
  write(u6,*) ' Sorry, ICISTR = 1 is out of fashion'
  write(u6,*) ' Switch to ICISTR = 2 - or reprogram'
  !stop ' DENSI2T : ICISTR = 1 in use'
  call SYSABENDMSG('lucia_util/densi2','Internal error','')
else if (ICISTR >= 2) then
  S2_TERM1 = Zero
  call GASDN2(I12,RHO1,RHO2,RHO2S,RHO2A,L,R,VEC3,NACOB,NSTSO(IATP)%A,NSTSO(IBTP)%A,NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA, &
              NOCTPB,NIRREP,NSMOB,MXPNGAS,NOBPTS,IOBPTS,MAXK,MAXI,VEC3(1+KCSCR),VEC3,NGAS,NELFSPGP,IDC,I1,XI1S,I2,XI2S,I3,XI3S,I4, &
              XI4S,INSCR,RHO1S,LUL,LUR,PSSIGN,PSSIGN,NBATCHL,LLBTL,LI1BTL,LIBTL,NBATCHR,LLBTR,LI1BTR,LIBTR,CONSPA,CONSPB,LSCLFCL, &
              LSCLFCR,S2_TERM1,IPHGAS,IDOSRHO1,SRHO1,IPACK)

  call GADSUM(RHO1,NACOB**2)
  if (I12 == 2) then
    if (IPACK) then
      ! If IPACK == .TRUE. then
      ! Number of elements in symmetric and antisymmetric 2-body
      ! density matrices are given in Nijkl.
      NIJ = nTri_Elem(NACOB)
      NIJKL = nTri_Elem(NIJ)
      call GADSUM(RHO2S,NIJKL)
      call GADSUM(RHO2A,NIJKL)
    else
      call GADSUM(RHO2,nTri_Elem(NACOB**2))
    end if
  end if
  if (IDOSRHO1 == 1) call GADSUM(SRHO1,NACOB**2)
  call GADSUM_SCAL(S2_TERM1)

  ! CALL GASDN2 --> 89
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
call NATORB_LUCIA(RHO1,NSMOB,NTOOBS,NACOBS,NINOBS,IREOST,XNATO,RHO1SM,OCCSM,NACOB,RHO1P)

if (IPRDEN >= 5) then
  write(u6,*) ' One-electron density matrix'
  write(u6,*) ' ==========================='
  call WRTMAT(RHO1,NTOOB,NTOOB,NTOOB,NTOOB)
  if (I12 == 2) then
    write(u6,*) ' Two-electron density'
    call PRSYM(RHO2,NACOB**2)
  end if
end if

if (I12 == 2) then
  ! <L!S**2|R>
  EXPS2 = S2_TERM1+Quart*real(4*NAEL+(NAEL-NBEL)*(NAEL-NBEL-2),kind=wp)
  if (IPRDEN > 0) then
    write(u6,*) ' Term 1 to S2 ',S2_TERM1
    write(u6,*) ' Expectation value of S2 ',EXPS2
  end if
else
  EXPS2 = Zero
end if

if ((IDOSRHO1 == 1) .and. (IPRDEN >= 2)) then
  write(u6,*) ' One-electron spindensity <0!E(aa) - E(bb)!0>'
  call WRTMAT(SRHO1,NTOOB,NTOOB,NTOOB,NTOOB)
end if

! Eliminate local memory
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
call mma_deallocate(OCSTR)
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

end subroutine DENSI2
