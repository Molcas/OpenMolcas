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
! Copyright (C) 1994-1996, Jeppe Olsen                                 *
!***********************************************************************

subroutine DENSI2_MCLR(I12,RHO1,RHO2,L,R,LUL,LUR,ieaw,n1,n2)
! Density matrices between L and R
!
! I12 = 1 => only one-body density
! I12 = 2 => one- and two-body density matrices
!
! Jeppe Olsen,      Oct 1994
! GAS modifications Aug 1995
! Two body density added, 1996
!
! Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
! ijkl = ij*(ij-1)/2+kl, ij >= kl
!
! =====
! Input
! =====
! Definition of L and R is picked up from CANDS
! with L being S and  R being C

use Index_Functions, only: nTri_Elem
use Str_Info, only: IATPM1, IATPM2, IBTPM1, IBTPM2, ITYP_DUMMY, NELEC, NOCTYP, STR
use MCLR_Data, only: IASTFI, IBSTFI, IBTSOB, ICISTR, IDC, MAXI, MAXK, MNR1IC, MXINKA, MXPNGAS, MXR3IC, MXSB, MXSOOB, NACOB, &
                     NOBPTS, NTSOB, PSSIGN
use CandS, only: ICSM, ICSPC, ISSM, ISSPC
use input_mclr, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: I12, LUL, LUR, ieaw, n1, n2
real(kind=wp), intent(out) :: RHO1(NACOB,NACOB), RHO2(nTri_Elem(NACOB**2))
real(kind=wp), intent(inout) :: L(*), R(*)
integer(kind=iwp) :: IATP, IATP2, IBTP, IBTP2, idum(1), IDUMMY(1), IMNMX, INTSCR, IOCTPA, IOCTPB, JATP, JBTP, LSCR1, LSCR12, &
                     LSCR2, LSCR3, MAXA, MAXA1, MAXB, MAXB1, MAXIK, MXCIJA, MXCIJAB, MXCIJB, MXCJ, MXIJST, MXIJSTF, MXSTBL, &
                     MXSXBL, MXTSOB, NAEL, NBEL, NGAS, NOCTPA, NOCTPB, NOOS
integer(kind=iwp), allocatable :: CBLTP(:), CIOIO(:), IX(:,:), OOS(:,:), SBLTP(:), SIOIO(:)
real(kind=wp), allocatable :: C2(:), CB(:), INSCR(:), RHO1P(:), RHO1S(:), SB(:), XIXS(:,:), XNATO(:)

IDUM = 0
NGAS = 3

RHO1(:,:) = Zero
RHO2(:) = Zero

! Info for this internal space

IATP = IASTFI(ISSPC)
IBTP = IBSTFI(ISSPC)
JATP = IASTFI(ICSPC)
JBTP = IBSTFI(ICSPC)
if ((IATP /= JATP) .or. (IBTP /= JBTP)) then
  write(u6,*) ' My world is falling apart'
  write(u6,*) ' C and sigma belongs to different types of strings'
  write(u6,*) ' IATP IBTP JATP JBTP ',IATP,IBTP,JATP,JBTP
  call Abend()
end if
NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
! Offsets for supergroups
IOCTPA = 1
IOCTPB = 1

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)

! Largest number of strings of given symmetry and type
MAXA = 0
if (NAEL >= 1) then
  MAXA1 = IMNMX(Str(IATPM1)%NSTSO,nIrrep*NOCTYP(IATPM1),2)
  MAXA = max(MAXA,MAXA1)
end if
if (NAEL >= 2) then
  MAXA1 = IMNMX(Str(IATPM2)%NSTSO,nIrrep*NOCTYP(IATPM2),2)
  MAXA = max(MAXA,MAXA1)
end if
MAXB = 0
if (NBEL >= 1) then
  MAXB1 = IMNMX(Str(IBTPM1)%NSTSO,nIrrep*NOCTYP(IBTPM1),2)
  MAXB = max(MAXB,MAXB1)
end if
if (NBEL >= 2) then
  MAXB1 = IMNMX(Str(IBTPM2)%NSTSO,nIrrep*NOCTYP(IBTPM2),2)
  MAXB = max(MAXB,MAXB1)
end if
MXSTBL = max(MAXA,MAXB)
! Largest number of resolution strings and spectator strings
! that can be treated simultaneously
! replace with MXINKA !!!
MAXI = min(MXINKA,MXSTBL)
MAXK = min(MXINKA,MXSTBL)
! Largest active orbital block belonging to given type and symmetry
MXTSOB = max(0,maxval(NOBPTS(1:NGAS,1:nIrrep)))
! Local scratch arrays for blocks of C and sigma
LSCR1 = 0
if (ICISTR <= 2) then
  LSCR1 = MXSB
else if (ICISTR == 3) then
  LSCR1 = MXSOOB
end if
if (ICISTR == 1) then
  call mma_allocate(CB,LSCR1,Label='CB')
  call mma_allocate(SB,LSCR1,Label='SB')
end if

! SCRATCH space for block of two-electron density matrix
! A 4 index block with four indices belonging OS class

INTSCR = MXTSOB**4

call mma_allocate(INSCR,INTSCR,Label='INSCR')

! Arrays giving allowed type combinations
call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')

call IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,SIOIO)

call IAIBCM_MCLR(MNR1IC(ICSPC),MXR3IC(ICSPC),NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,CIOIO)

! Get memory requirements

IATP2 = min(IATP+2,ITYP_Dummy)
IBTP2 = min(IBTP+2,ITYP_Dummy)
call MXRESC(CIOIO,IATP,IBTP,NOCTPA,NOCTPB,nIrrep,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,Str(IATP+1)%NSTSO,NOCTYP(IATP+1), &
            Str(IBTP+1)%NSTSO,NOCTYP(IBTP+1),NTSOB,MAXK,Str(IATP2)%NSTSO,NOCTYP(IATP2),Str(IBTP2)%NSTSO,NOCTYP(IBTP2), &
            Str(IATP)%EL123,Str(IBTP)%EL123,MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,MXIJSTF)

LSCR2 = max(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
LSCR12 = max(LSCR1,2*LSCR2)
call mma_allocate(C2,LSCR12,Label='C2')

! Space for annihilation/creation mappings
MAXIK = max(MAXI,MAXK)
LSCR3 = max(MXSTBL*MXTSOB,MXIJST,MAXIK*MXTSOB*MXTSOB)
call mma_allocate(IX,LSCR3,4,Label='IX')
call mma_allocate(XIXS,LSCR3,4,Label='XIXS')
! Arrays giving block type
call mma_allocate(SBLTP,nIrrep,Label='SBLTP')
call mma_allocate(CBLTP,nIrrep,Label='CBLTP')
! Arrays for additional symmetry operation
call ZBLTP(ISSM,nIrrep,IDC,SBLTP,idum)
call ZBLTP(ICSM,nIrrep,IDC,CBLTP,idum)
!10 OOS array
NOOS = NOCTPA*NOCTPB*nIrrep
call mma_allocate(OOS,NOOS,10,Label='OSS')
! scratch space containing active one body
call mma_allocate(RHO1S,NACOB**2,Label='RHO1S')
! For natural orbitals
call mma_allocate(RHO1P,nTri_Elem(NACOB),Label='RHO1P')
call mma_allocate(XNATO,NACOB**2,Label='XNATO')
! Natural orbitals in symmetry blocks
!call mma_allocate(RHO1SM,NACOB ** 2,Label='RHO1SM')
!call mma_allocate(XNATSM,NACOB ** 2,Label='XNATSM')
!call mma_allocate(OCCSM,NACOB,Label='OCCSM')

! Transform from combination scaling to determinant scaling

if ((IDC /= 1) .and. (ICISTR == 1)) then
  ! Left CI vector
  call SCDTC2_MCLR(L,ISSM,SBLTP,nIrrep,NOCTPA,NOCTPB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,SIOIO,IDC,2,IDUMMY)
  ! Right CI vector
  call SCDTC2_MCLR(R,ICSM,CBLTP,nIrrep,NOCTPA,NOCTPB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,CIOIO,IDC,2,IDUMMY)
end if

if (ICISTR == 1) then
  call GASDN2_MCLR(I12,RHO1,RHO2,R,L,CB,SB,C2,CIOIO,SIOIO,ICSM,ISSM,CBLTP,SBLTP,NACOB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,NAEL,IATP, &
                   NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,nIrrep,MXPNGAS,NTSOB,IBTSOB,MAXK,MAXI,LSCR1,LSCR1,C2(LSCR2+1:LSCR12), &
                   C2(1:LSCR2),NGAS,Str(IATP)%EL123,Str(IBTP)%EL123,IDC,OOS(:,1),OOS(:,2),OOS(:,3),OOS(:,4),OOS(:,5),OOS(:,6), &
                   OOS(:,7),OOS(:,8),OOS(:,9),OOS(:,10),IX(:,1),XIXS(:,1),IX(:,2),XIXS(:,2),IX(:,3),XIXS(:,3),IX(:,4),XIXS(:,4), &
                   INSCR,RHO1S,LUL,LUR,PSSIGN,PSSIGN,ieaw,n1,n2)
else if (ICISTR >= 2) then
  call GASDN2_MCLR(I12,RHO1,RHO2,R,L,R,L,C2,CIOIO,SIOIO,ICSM,ISSM,CBLTP,SBLTP,NACOB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,NAEL,IATP, &
                   NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,nIrrep,MXPNGAS,NTSOB,IBTSOB,MAXK,MAXI,LSCR1,LSCR1,C2(LSCR2+1:LSCR12), &
                   C2(1:LSCR2),NGAS,Str(IATP)%EL123,Str(IBTP)%EL123,IDC,OOS(:,1),OOS(:,2),OOS(:,3),OOS(:,4),OOS(:,5),OOS(:,6), &
                   OOS(:,7),OOS(:,8),OOS(:,9),OOS(:,10),IX(:,1),XIXS(:,1),IX(:,2),XIXS(:,2),IX(:,3),XIXS(:,3),IX(:,4),XIXS(:,4), &
                   INSCR,RHO1S,LUL,LUR,PSSIGN,PSSIGN,ieaw,n1,n2)
end if

if ((IDC /= 1) .and. (ICISTR == 1)) then
  ! Transform from combination scaling to determinant scaling

  call SCDTC2_MCLR(L,ISSM,SBLTP,nIrrep,NOCTPA,NOCTPB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,SIOIO,IDC,1,IDUMMY)
  call SCDTC2_MCLR(R,ICSM,CBLTP,nIrrep,NOCTPA,NOCTPB,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,CIOIO,IDC,1,IDUMMY)
end if

! Free memory

if (ICISTR == 1) then
  call mma_deallocate(CB)
  call mma_deallocate(SB)
end if
call mma_deallocate(INSCR)
call mma_deallocate(SIOIO)
call mma_deallocate(CIOIO)
call mma_deallocate(C2)
call mma_deallocate(IX)
call mma_deallocate(XIXS)
call mma_deallocate(SBLTP)
call mma_deallocate(CBLTP)
call mma_deallocate(OOS)
call mma_deallocate(RHO1S)
call mma_deallocate(RHO1P)
call mma_deallocate(XNATO)

end subroutine DENSI2_MCLR
