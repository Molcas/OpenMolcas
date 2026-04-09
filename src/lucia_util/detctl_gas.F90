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

subroutine DETCTL_GAS()

use CandS, only: ICSM, ICSPC, ISSM, ISSPC
use lucia_data, only: Allocate_Local_Arrays, CBLTP, CI1BT, CIBT, CLBT, CLEBT, CONF_OCC, Deallocate_Local_Arrays, IADVICE, &
                      IBSPGPFTP, IDC, IGSOCC, IPHGAS, IPRCIX, IREFSM, kvec3_length, LCSBLK, MAX_STR_OC_BLK, MAX_STR_SPGP, &
                      Memory_Needed_Lucia, MNHL, MXINKA, MXNSTR, MXNTTS, MXPCSM, MXPNGAS, MXPNSMST, MXPORB, MXSOOB, MXTSOB, &
                      NCONF_PER_OPEN, NCONF_PER_SYM, NCSF_HEXS, NCSF_PER_SYM, NELEC, NELFSPGP, NELFTP, NGAS, NHLFSPGP, NIRREP, &
                      NOBPT, NOBPTS, NOCOB, NOCSF, NOCTYP, NPCSCNF, NPDTCNF, NSD_PER_SYM, NSMOB, NSTFSMSPGP, NSTSO, NTOOB, PSSIGN, &
                      SDREO_I, XISPSM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: IATP, IATPM1, IATPM2, IBASSPC(1), IBTP, IBTPM1, IBTPM2, INTSCR, IOCCLS(1), IOCTPA, IOCTPB, JSYM, LBLOCK, &
                     LSCR12, LSCR2, LSCR3, LZ, LZSCR, MAXA, MAXA1, MAXB, MAXB1, MAXI, MAXIK, MAXK, MX_NSPII, MXADKBLK, &
                     MXADKBLK_AS, MXCIJA, MXCIJAB, MXCIJB, MXCJ, MXCJ_ALLSYM, MXSTBL, MXSTBL0, MXSXBL, NAEL, NBATCH, NBEL, NBLOCK, &
                     NDET, NEL, NOCCLS, NOCTPA, NOCTPB, NTTS
integer(kind=iwp), allocatable :: BASSPC(:), KLOCCLS(:), LCIOIO(:), SVST(:)
integer(kind=iwp), external :: IMNMX

! Set variables in Module cands
JSYM = IREFSM
ICSM = JSYM
ISSM = JSYM
ICSPC = 1
ISSPC = 1
! Set NDET
NDET = int(XISPSM(JSYM,1))
! Number of occupation classes
IATP = 1
IBTP = 2
NEL = NELFTP(IATP)+NELFTP(IBTP)
call OCCLS(1,NOCCLS,IOCCLS,NEL,NGAS,IGSOCC(:,1),IGSOCC(:,2),0,IBASSPC,NOBPT)
! and then the occupation classes
call mma_allocate(KLOCCLS,NGAS*NOCCLS,Label='KLOCCLS')
call mma_allocate(BASSPC,NOCCLS,Label='BASSPC')
call OCCLS(2,NOCCLS,KLOCCLS,NEL,NGAS,IGSOCC(:,1),IGSOCC(:,2),1,BASSPC,NOBPT)
call mma_deallocate(BASSPC)
if (NOCSF == 0) then
  ! Initial information on CSF expansion
  call CSFDIM_GAS(KLOCCLS,NOCCLS,JSYM)
  ! Prototype dets and csf's and CSF'SD matrices
  call CSDTMT_GAS()
end if
! Allocate memory for diagonalization
!if (ISIMSYM == 0) then
LBLOCK = MXSOOB
!else
!  LBLOCK = MXSOOB_AS
!end if
LBLOCK = max(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
!PAM06 LBLOCK = MAX(XISPSM(IREFSM,1),real(MXSOOB,kind=wp))
LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))

! Information about block structure- needed by new PICO2 routine.
! Memory for partitioning of C vector
IATP = 1
IBTP = 2
NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
NTTS = MXNTTS
call Allocate_Local_Arrays(NTTS,NIRREP)
! Additional info required to construct partitioning

call mma_allocate(LCIOIO,NOCTPA*NOCTPB,Label='LCIOIO')

call IAIBCM(ICSPC,LCIOIO)
call mma_allocate(SVST,1,Label='SVST')
call ZBLTP(jsym,NIRREP,IDC,CBLTP,SVST)
call mma_deallocate(SVST)

! Batches  of C vector
call PART_CIV2(IDC,NSTSO(IATP)%A,NSTSO(IBTP)%A,NOCTPA,NOCTPB,NIRREP,LCIOIO,jsym,NBATCH,CLBT,CLEBT,CI1BT,CIBT,0)
! Number of BLOCKS
NBLOCK = CI1BT(NBATCH)+CLBT(NBATCH)-1
! Length of each block
call EXTRROW(CIBT,8,8,NBLOCK,CI1BT)
if (NEL > 0) call CNFORD_GAS(KLOCCLS,NOCCLS,jsym,SDREO_I(jsym)%A,CIBT,NBLOCK)

call Deallocate_Local_Arrays()
! If PICO2/SBLOCK are used, three blocks are used in PICO2, so
!
! Sblock is used in general nowadays so,
! Largest block of strings in zero order space
MXSTBL0 = MXNSTR
! type of alpha and beta strings
IATP = 1
IBTP = 2
! alpha and beta strings with an electron removed
IATPM1 = 3
IBTPM1 = 4
! alpha and beta strings with two electrons removed
IATPM2 = 5
IBTPM2 = 6

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
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

MXSTBL = max(MAXA,MAXB,MXSTBL0)

if (IPRCIX >= 2) write(u6,*) ' Largest block of strings with given symmetry and type',MXSTBL
! Largest number of resolution strings and spectator strings
! that can be treated simultaneously
MAXI = min(MXINKA,MXSTBL)
MAXK = min(MXINKA,MXSTBL)
! scratch space for projected matrices and a CI block

! Scratch space for CJKAIB resolution matrices
! Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)
call MXRESCPH(LCIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NIRREP,NSTFSMSPGP,MXPNSMST,NSMOB,MXPNGAS,NGAS,NOBPTS,MAXK,NELFSPGP,MXCJ,MXCIJA, &
              MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)

call mma_deallocate(LCIOIO)

if (IPRCIX >= 2) then
  write(u6,*) 'DETCTL : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
  write(u6,*) ' MXADKBLK, MXADKBLK_AS',MXADKBLK,MXADKBLK_AS
end if
!if (ISIMSYM == 1) then
!  MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
!  MXADKBLK_AS = MXADKBLK
!end if
! Using hardwired routines, MXCIJAB also used
LSCR2 = max(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
if (IPRCIX >= 2) write(u6,*) ' Space for two resolution matrices ',2*LSCR2
LSCR12 = max(LBLOCK,2*LSCR2)

KVEC3_LENGTH = max(LSCR12,2*LBLOCK)

! Calculate how much memory the sigma routines needs.

! Memory needed in Sigma_Master
MEMORY_NEEDED_LUCIA = KVEC3_LENGTH
! Memory needed in MV7
MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA+2*1+NOCTPA*NOCTPB+NIRREP+11*MXNTTS
!if ((IDC == 3) .or. (IDC == 4)) MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA+NIRREP
! Memory needed in SBLOCK
! Largest active orbital block belonging to given type and symmetry
MXTSOB = max(0,maxval(NOBPTS(1:NGAS,1:NSMOB)))
INTSCR = max(MXTSOB**4,NTOOB**2)
! vectors able to hold strings of given sym and type
MAXIK = max(MAXI,MAXK)
LSCR3 = max(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXNSTR)
! Space for four blocks of string occupations and arrays of
! reordering arrays
! Also used to hold an NORB*NORB matrix
LZSCR = (max(NAEL,NBEL)+3)*(NOCOB+1)+2*NOCOB+NOCOB*NOCOB
LZ = (max(NAEL,NBEL)+2)*NOCOB
MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA+NOCTPA**2+NOCTPB**2+2*NIRREP**2+2*INTSCR+3*NOCTPA*NOCTPB+NIRREP+8*LSCR3+12*NTTS+ &
                      MAX_STR_OC_BLK+2*MAX_STR_SPGP+2*LZ+LZSCR
!if ((IDC == 3) .or. (IDC == 4))  MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA+NIRREP
!if (ISIMSYM == 0) then
MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA+4*MAX_STR_SPGP
!end iF

call LUCIA2MOLCAS(CONF_OCC(jsym)%A,SDREO_I(jsym)%A,ndet,ncsf_per_sym,nsd_per_sym,nconf_per_sym,mxpcsm,mxporb,nconf_per_open, &
                  npdtcnf,npcscnf,nCSF_HEXS)

call mma_deallocate(KLOCCLS)

end subroutine DETCTL_GAS
