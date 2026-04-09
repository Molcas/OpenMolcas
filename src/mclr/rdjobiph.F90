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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUGPRINT_
subroutine RdJobIph(CIVec)
!***********************************************************************
!                                                                      *
!     Read the contents of the JOBIPH file.                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use MckDat, only: sNew
use gugx, only: CIS, EXS, SGS
use MCLR_Data, only: CMO, FnJob, FnMck, G1t, G2t, IRLXROOT, ISNAC, ISTATE, LuJob, LuMck, nA, NACSTATES, nNA, NSSA, OVERRIDE, SA
use input_mclr, only: Debug, ERASSCF, Headerjp, iMCPD, iMSPD, iPT2, iRoot, iSpin, iTOC, iTocIph, lRoots, McKinley, nActEl, nAsh, &
                      nBas, nCOnf, nCSF, nDel, nElec3, nFro, nHole1, nIsh, nOrb, nRoots, nRS1, nRS2, nRS3, nSym, ntAsh, ntASqr, &
                      ntATri, ntBas, ntBSqr, ntBTri, ntIsh, ntISqr, ntITri, PT2, State_Sym, TitleJP, Weight
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use Molcas, only: LenIn, MxOrb, MxRoot, MxSym
use RASDim, only: MxIter, MxTit
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), allocatable, intent(out) :: CIVec(:,:)
integer(kind=iwp) :: i, iDisk, iDummer, iGo, iMode, iOpt, iRC, iS, iSym, Iter, j, jS, k, kRoots, kS, Length, lS, nAct, nAct2, &
                     nAct4, nG1, nG2
real(kind=wp) :: PotNuc0, rdum(1), Temp
logical(kind=iwp) :: Found
character(len=8) :: Method
real(kind=wp), allocatable :: Tmp2(:)
character, allocatable :: TempTxt(:)

!                                                                      *
!***********************************************************************
!                                                                      *
debug = .false.
#ifdef _DEBUGPRINT_
debug = .true.
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
!     Save the ROOT input parameter                                    *
!----------------------------------------------------------------------*
kRoots = lRoots
!----------------------------------------------------------------------*
!     Read the table of disk addresses                                 *
!----------------------------------------------------------------------*
call DaName(LuJob,FnJob)
iDisk = 0
call iDaFile(LuJob,2,iToc,iTOCIPH,iDisk)
!----------------------------------------------------------------------*
!     Read the the system description                                  *
!----------------------------------------------------------------------*
call mma_allocate(TempTxt,(LenIn+8)*MxOrb,Label='TempTxt')
iDisk = iToc(1)

!write(u6,*) 'if dmrg, it should be something else'
call WR_RASSCF_Info(LuJob,2,iDisk,nActEl,iSpin,nSym,State_sym,nFro,nIsh,nAsh,nDel,nBas,MxSym,TempTxt,(LenIn+8)*mxorb,nConf, &
                    HeaderJP,144,TitleJP,4*18*mxTit,PotNuc0,lRoots,nRoots,iRoot,mxRoot,nRs1,nRs2,nRs3,nHole1,nElec3,iPt2,Weight)

if (doDMRG) then ! yma
  nash(:) = LRras2(:)
  nrs2(:) = LRras2(:)
end if
!do i=1,8
!  write(u6,*) i,'-irrep',nIsh(i),nAsh(i),nRs1(i),nRs2(i),nRs3(i)
!  call xflush(u6)
!end do

call mma_deallocate(TempTxt)
!----------------------------------------------------------------------*
!     Overwrite the variable lroots if approriate, i.e if lroot        *
!     was set by input.                                                *
!----------------------------------------------------------------------*
if (kRoots /= -1) then
  if (iPt2 /= 0) then
    write(u6,*) 'RdJobiph: kRoots /= -1 .and. iPt2 /= 0'
    call Abend()
  else if (kRoots > lRoots) then
    write(u6,*) 'RdJobiph: kRoots /= -1 .and. kRoots > lRoots'
    call Abend()
  end if
  lRoots = kRoots
  nRoots = 1
end if
!----------------------------------------------------------------------*
!     Precompute the total sum of variables and size of matrices       *
!----------------------------------------------------------------------*
ntIsh = sum(nIsh(1:nSym))
ntIsqr = sum(nIsh(1:nSym)**2)
ntAsh = sum(nAsh(1:nSym))
ntAsqr = sum(nAsh(1:nSym)**2)
ntBas = sum(nBas(1:nSym))
ntBsqr = sum(nBas(1:nSym)**2)
nOrb(1:nSym) = nBas(1:nSym)-nDel(1:nSym)
Length = sum(nBas(1:nSym)*nOrb(1:nSym))
ntItri = 0
ntAtri = 0
ntBtri = 0
nna = 0
do iSym=1,nSym
  ntItri = ntItri+nTri_Elem(nIsh(iSym))
  ntAtri = ntAtri+nTri_Elem(nAsh(iSym))
  ntBtri = ntBtri+nTri_Elem(nBas(iSym))
  nA(iSym) = nna
  nnA = nnA+nAsh(isym)
end do

! Generate the nr. of csf in each sub-sym, used in geom-opt with SA DMRG-SCF
if (doDMRG) then  ! yma
  imode = -99
  ! generate the Nr. of csfs in each sym
  call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,rdum,imode,State_Sym,State_Sym)
  NCSF(1:nSym) = CIS%NCSF(1:nSym)
  NCONF = CIS%NCSF(State_Sym)
  call mkGuga_Free(SGS,CIS,EXS)

  !do isym=1,8
  !  write(u6,*) 'isym_ncsf in rdjobiph ',ncsf(isym)
  !end do
end if

!----------------------------------------------------------------------*
!     Load the orbitals used in the last macro iteration               *
!----------------------------------------------------------------------*

call mma_allocate(CMO,Length,Label='CMO')
call Get_dArray_chk('Last orbitals',CMO,Length)

! Read state for geo opt

call Get_iScalar('Relax CASSCF root',irlxroot)
call Get_cArray('Relax Method',Method,8)
iMCPD = .false.
iMSPD = .false.
if ((Method == 'MCPDFT') .or. (Method == 'MSPDFT')) then
  iMCPD = .true.
  if (Method == 'MSPDFT') iMSPD = .true.
  do i=1,lroots
    if (iroot(i) == irlxroot) istate = i
  end do
end if
if ((Method == 'CASSCFSA') .or. (Method == 'CASPT2') .or. (Method == 'RASSCFSA')) then
  call Get_iScalar('SA ready',iGo)
  if (iGo /= 2) SA = .true.
  Found = .true.
  if (override) then
    if (isNAC) then
      do j=1,2
        NSSA(j) = 0
        do i=1,lroots
          if (iroot(i) == NACStates(j)) NSSA(j) = i
        end do
        if (NSSA(j) == 0) Found = .false.
      end do
    else
      irlxroot = iroot(istate)
    end if
  else
    istate = 0
    do i=1,lroots
      if (iroot(i) == irlxroot) istate = i
    end do
    if (istate == 0) Found = .false.
  end if
  if (.not. Found) then
    call WarningMessage(2,'Cannot relax a root not included in the SA')
    call Abend()
  end if
else if ((irlxroot == 1) .and. (.not. (McKinley .or. PT2 .or. iMCPD))) then
  write(u6,*)
  write(u6,*) 'W A R N I N G !'
  write(u6,*)
  write(u6,*) 'Redundant rlxroot input in RASSCF!'
  write(u6,*) 'I''ll sign off here without a clean termination!'
  write(u6,*) 'However, I have to fix the epilogue file.'
  write(u6,*)
  irc = -1
  iopt = ibset(0,sNew)
  call OPNMCK(irc,iopt,FNMCK,LUMCK)
  iopt = 0
  call WrMck(iRC,iOpt,'nSym',1,nBas,iDummer)
  call ClsFls_MCLR()
  call Finish(0)
end if

!iDisk = iToc(9)
!if (IPT2 == 0) iDisk = iToc(2)
!call dDaFile(LuJob,2,CMO,ntBsqr,iDisk)
!jpCMO = 1
!do iSym=1,nSym
!  CMO(jpCMO+norb(isym)*nbas(isym):jpCMO+(norb(isym)+ndel(isym))*nbas(isym)-1) = Zero
!  write(Line,'(A,i2.2)') 'MO coefficients, iSym = ',iSym
!  call RecPrt(Line,' ',CMO(jpCMO),nBas(iSym),nBas(iSym))
!  jpCMO = jpCMO+nBas(iSym)*nBas(iSym)
!end do
!----------------------------------------------------------------------*
!     Load the CI vectors for the SA roots                             *
!----------------------------------------------------------------------*

! If doDMRG, introducing CI coeffieients later :
!    1) only coefficients of importants DETs using MPS2CI
!    2) and together with DET numbers from GUGA generation part

if (doDMRG) then  ! yma
  call mma_allocate(CIVec,nConf,nroots,Label='CIVec')
else
  call mma_allocate(CIVec,nConf,nroots,Label='CIVec')
  do i=1,nroots
    j = iroot(i)
    iDisk = iToc(4)
    do k=1,j-1
      call dDaFile(LuJob,0,rdum,nConf,iDisk)
    end do
    call dDaFile(LuJob,2,CIVec(:,i),nConf,iDisk)
  end do
!#ifdef _DEBUGPRINT_       ! yma umcomment
!  do i=1,nroots            !yma
!    inum = 0
!    dv_ci2 = Zero
!    do j=1,nconf
!      !yma CI-threshold
!      if (abs(CIVec(j,i)) < Zero) then
!        inum = inum+1
!        CIVec(j,i) = Zero
!      else
!        dv_ci2 = dv_ci2+CIVec(j,i)**2
!      end if
!    end do
!    !call DVcPrt('CI coefficients',' ',CIVec(:,i),nConf)!yma
!    !write(u6,*) 'dismissed dets num', inum
!    !write(u6,*) 'absolutely CI^2',dv_ci2
!  end do
!#endif
end if

!----------------------------------------------------------------------*
!     Load state energy                                                *
!----------------------------------------------------------------------*
call mma_allocate(Tmp2,mxRoot*mxIter,Label='Tmp2')
iDisk = iToc(6)
#ifdef _DEBUGPRINT_
if (debug) then
  write(u6,*) 'NROOTS: ',nroots
  write(u6,*) 'iROOTS: ',(iroot(i),i=1,nroots)
  write(u6,*) 'lROOTS: ',lroots
end if
#endif
call dDaFile(LuJob,2,Tmp2,mxRoot*mxIter,iDisk)

do iter=0,mxIter-1
  do i=1,nroots
    j = iroot(i)
    ! It should be 0.0 in DMRG case
    Temp = Tmp2(iter*mxRoot+j)
    if (Temp /= Zero) ERASSCF(i) = Temp
    !if (debug) write(u6,*) ERASSCF(i),i
  end do
end do

#ifdef _DEBUGPRINT_
if (debug) then
  write(u6,*) (Tmp2(i),i=1,lroots)
  write(u6,*) 'RASSCF energies=',(ERASSCF(i),i=1,nroots)
end if
#endif
call mma_deallocate(Tmp2)

nAct = 0    ! 1/2

if (doDMRG) then  ! yma
  call dmrg_dim_change_mclr(RGras2,ntash,0)
  nash(:) = RGras2(:)
  nrs2(:) = RGras2(:)
end if

nAct2 = 0
nAct4 = 0
do iSym=1,nSym
  nAct = nAct+nAsh(iSym)
  nAct2 = nAct2+nAsh(iSym)**2
end do
do iS=1,nSym
  do jS=1,nSym
    do kS=1,nSym
      lS = Mul(Mul(is,js),ks)
      nAct4 = nAct4+nAsh(iS)*nAsh(jS)*nAsh(kS)*nAsh(lS)
    end do
  end do
end do

nG1 = nTri_Elem(nAct)
call mma_allocate(G1t,nG1,Label='G1t')
nG2 = nTri_Elem(nG1)
call mma_allocate(G2t,nG2,Label='G2t')
call RDDENS(G1t,ng1,G2t,ng2)

#ifdef _DEBUGPRINT_
call Triprt('G1',' ',G1t,ntash)
call Triprt('G2',' ',G2t,ng1)
#endif

if (doDMRG) then ! yma
  call dmrg_dim_change_mclr(LRras2,nact,0)
  nash(:) = LRras2(:)
  nrs2(:) = LRras2(:)
end if
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine RdJobIph
