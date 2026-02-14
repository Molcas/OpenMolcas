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
! Copyright (C) 1989, Bjorn O. Roos                                    *
!               1989, Per Ake Malmqvist                                *
!               1991,1993,1996, Markus P. Fuelscher                    *
!***********************************************************************

subroutine SXCtl(CMO,OCC,D,P,PA,FI,FA,D1A,THMAX,IFINAL)
!*******************************************************************
!                                                                  *
! super-CI control section                                         *
!                                                                  *
! calling arguments:                                               *
! CMO     : array of real*8                                        *
!           MO-coefficients                                        *
! OCC     : array of real*8                                        *
!           orbital occupation numbers                             *
! D       : array of real*8                                        *
!           averaged one-body density matrix                       *
! P       : array of real*8                                        *
!           averaged two-body density matrix                       *
! PA      : array of real*8                                        *
!           averaged antisymmetric two-body density matrix         *
! FI      : array of real*8                                        *
!           inactive Fock matrix                                   *
! FA      : array of real*8                                        *
!           active Fock matrix                                     *
! D1A     : array of real*8                                        *
!           active one-body density matrix in AO basis             *
! IFINAL  : integer                                                *
!           termination flag                                       *
!                                                                  *
!------------------------------------------------------------------*
!                                                                  *
! written by:                                                      *
! B.O. Roos and P.Aa. Malmqvist                                    *
! University of Lund, Sweden, 1989                                 *
!                                                                  *
!------------------------------------------------------------------*
!                                                                  *
! history:                                                         *
! - updated for MOLCAS version 2                                   *
!   M.P. Fuelscher, University of Lund, Sweden, 1991               *
! - updated for MOLCAS version 3                                   *
!   M.P. Fuelscher, University of Lund, Sweden, 1993               *
! - updated for integral direct and reaction field calculations    *
!   M.P. Fuelscher, University of Lund, Sweden, 1996               *
!                                                                  *
!*******************************************************************

use fciqmc, only: DoNECI
use Fock_util_global, only: ALGO, DoCholesky
use Lucia_Interface, only: Lucia_Util
use wadr, only: BM, DIA, F1, F2, NLX, SXG, SXH, SXN
use input_ras, only: KeyHEUR
use rasscf_global, only: DoBlockDMRG, DoDMRG, ECAS, EMY, ESX, ExFac, IADR15, iCIOnly, iPT2, ISTORP, ITER, ITERSX, ITMAX, KSDFT, &
                         l_casdft, NAC, nDimSX, nFint, NO2M, nQune, NROOT, NSXS, NTOT4, QNSTEP, QNUPDT, SXSEL, TMIN, VIA
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: JOBIPH, LUINTM, LUQUNE, NACTEL, NASH, NBAS, NDEL, NFRO, NISH, NORB, NRS1, NRS2, NRS3, NSSH, NSYM, NTOT, &
                        NTOT1, NTOT2
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use raswfn, only: wfn_mocoef, wfn_occnum
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO(*), OCC(*), D(*), P(*), PA(*), FI(*), FA(*), D1A(*), THMAX
integer(kind=iwp) :: IFINAL
integer(kind=iwp) :: i, iBas, IC, iDisk, IndT, IndType(56), iOff, iOrb, iPrLev, iRC, iRef, iShift, iSym, iWay, j, kMax, LCSXI, &
                     LuvvVec, MIAAE, MNO, NA, NAE, NAEAE, NAOAE, NB, NCR, NCR1, NE, NI, NIA, NIAIA, NLCC, NLHH, NLOVL, NLQ, NLX1, &
                     NLX2, NO, nP2Act, NQ
real(kind=wp) :: CASDFT_En, CIDUMMY(1), CPES, CPTS, Dummy(1), P2act(1), P2reo_size, TIOES, TIOS, XSXMAX
character(len=80) :: VecTyp
logical(kind=iwp) :: TraOnly
integer(kind=iwp), save :: nCall
real(kind=wp), allocatable :: CC(:), CMON(:), CMOX(:), CSX(:), DA(:), EDum(:), ENER_X(:), Fck(:), FTR(:), HH(:), OVL(:), P2Raw(:), &
                              P2reo(:), PUVX(:), QMat(:), QQ(:), SC(:), SCR(:), Sigma(:), SQ(:), STRP(:), SXDD(:), SXDF(:), &
                              SXHD(:), V1(:), V2(:), Vec(:), VL(:), VT(:), WO(:), X2(:), XMAT(:), XQN(:)
real(kind=wp), allocatable, target :: SMAT(:)
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: DDot_

! PAM01 The SXCI part has been slightly modified by P-AA M Jan 15, 2001:
! Changes affect several of the subroutines of this part.
! The program should now behave more gracefully when there are (almost)
! redundant orbital rotation modes.
! For individual 2x2 rotations that are to be excluded, the SX hamiltonian
! diagonal, used in the preconditioning step in DAVCRE, is set to a huge
! value. Corresponding elements of trial vectors and sigma vectors are
! zeroed.
! To take care of 'hidden' redundant linear combinations, a small extra
! term, proportional to the square norm of the orbital rotations, is
! added to the optimization problem: this means in practice that a small
! extra term is added in the SIGVEC routine.
! Also, an extra term is added to the overlap calculation. The
! proportionality factor in SIGVEC should be small, but larger than the
! one used in the COVLP routine.
! Presently I try the following values:
! Huge elements in diagonal values: 1.0e32
! (Tested against 1.0e20 in IF-statements)
! Extra term in overlaps (COVLP, SXHAM): 1.0e-14
! Extra term in SIGVEC:                  1.0e-12

! Local print level (if any)
IPRLEV = IPRLOC(4)
!write(u6,*) 'Entering SXCTL!'
if (IPRLEV >= DEBUG) write(u6,*) ' Entering SXCTL'

! --- Check for Cholesky ---------------

!call DecideOnCholesky(DoCholesky)
! --------------------------------------

! compute constants needed for addressing
NLX1 = 0
NLX2 = 0
NQ = 0
NSXS = 0
NIAIA = 0
NAEAE = 0
NAOAE = 0
MNO = 0
MIAAE = 0
NLX = 0
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  NA = NASH(ISYM)
  NI = NISH(ISYM)
  NE = NSSH(ISYM)
  NO = NORB(ISYM)
  NIA = NI+NA
  NAE = NA+NE
  NSXS = NSXS+NIA*NAE
  NIAIA = NIAIA+NIA**2
  NAEAE = NAEAE+NAE**2
  NAOAE = NAOAE+NA*NAE
  NLX1 = max(NLX1,NO**2)
  NLX2 = max(NLX2,NB*NO)
  NQ = max(NQ,NA*NO)
  MNO = max(MNO,NO)
  MIAAE = max(MIAAE,NIA*NAE)
  NLX = max(NLX,NIA*NA,NIA*NAE)
end do
if (NQ < NIAIA) NQ = NIAIA
NROOT = 1
NDIMSX = NROOT+NSXS
!***********************************************************************
! load back two-electron integrals (pu|vx)
!***********************************************************************
if ((.not. DoCholesky) .or. (ALGO == 1)) then
  if (nFint > 0) then
    iDisk = 0
    call mma_allocate(PUVX,nFint,Label='PUVX')
    call DDaFile(LUINTM,2,PUVX,nFint,iDisk)
  else
    call mma_allocate(PUVX,1,Label='PUVX')
  end if
end if
if (IPRLEV >= DEBUG) then
  write(u6,*) 'PUVX integrals in SXCTL'
  call wrtmat(PUVX,1,nFInt,1,nFInt)
end if
!********************************************************************************
! update and transform the Fock matrices FI and FA in MO basis ----> Fmat routine
!********************************************************************************
if ((.not. DoCholesky) .or. (ALGO == 1)) then
  call Fmat(CMO,PUVX,D,D1A,FI,FA)

else if (ALGO == 2) then

  ! Inactive-active contribution to ECAS
  call mma_allocate(DA,nTot1,Label='DA')
  call Fold(nSym,nBas,D1A,DA) !get the packed DA
  VIA = dDot_(nTot1,FI,1,DA,1)
  ECAS = EMY+VIA
  if (IPRLEV >= DEBUG) then
    write(u6,'(A,ES20.10)') ' Total core energy:            ',EMY
    write(u6,'(A,ES20.10)') ' inactive-active interaction:  ',VIA
    write(u6,'(A,ES20.10)') ' CAS energy (core+interaction):',ECAS
  end if
  call mma_deallocate(DA)

  TraOnly = .true.
  call CHO_CAS_DRV(irc,CMO,D,FI,D1A,FA,Dummy,TraOnly)

  if (irc /= 0) then
    write(u6,*) 'SXCTL: Cho_cas_drv non-zero return code! rc= ',irc
    call abend()
  end if

else

  write(u6,*) 'SXCTL: Illegal Cholesky parameter ALGO= ',ALGO
  call abend()

end if
!***********************************************************************
! reorder the two-body density matrix P
!***********************************************************************
if (.not. l_casdft) then
  ! ISTORP(NSYM+1) represents the size of the 2-body density matrix,d(vwxy), with vwxy all active.
  ! the size is computed as NAP*NAQ*NRS (sum over all symmetries). If Sym_R = Sym_S then triangular
  ! form over NRS... with R >= S, rectanguar otherwise.
  if (ISTORP(NSYM+1) > 0) then
    call mma_allocate(STRP,ISTORP(NSYM+1),Label='STRP')

    call PMAT_RASSCF(P,STRP)

    if ((ExFac /= One) .and. (KSDFT(1:3) /= 'SCF')) then
      call mma_allocate(P2reo,ISTORP(NSYM+1),Label='P2reo')
      call Get_Temp('nP2Act  ',P2Act,1)
      nP2Act = int(P2Act(1))
      call mma_allocate(P2RAW,nP2Act,Label='P2Raw')
      call Get_Temp('P2_RAW  ',P2RAW,nP2Act)
      call PMAT_RASSCF(P2RAW,P2reo)
      call mma_deallocate(P2Raw)
      P2reo_size = real(ISTORP(NSYM+1),kind=wp)
      call Put_Temp('nP2reo  ',[P2reo_size],1)
      call Put_Temp('P2_reo  ',P2reo,ISTORP(NSYM+1))
      call mma_deallocate(P2reo)
    end if
  else
    call mma_allocate(STRP,1,Label='STRP')
  end if
else ! GLM-CASDFT
  ! ISTORP(NSYM+1) here represents the size of the Dvw*Dxy array (product of one-body
  ! density matrix,d(vwxy), with vwxy all active. The size is computed as NAP*NAQ*NRS
  ! (sum over all symmetries). If Sym_R = Sym_S then triangular form over NRS...
  ! with R >= S, rectanguar otherwise. Basically we will use same simmetry as for dvwxy.
  if (ISTORP(NSYM+1) > 0) then
    !write(u6,*)
    !write(u6,*) ' ---------------------'
    call mma_allocate(STRP,ISTORP(NSYM+1),Label='STRP')
    call DmatDmat(D,STRP)
  else
    call mma_allocate(STRP,1,Label='STRP')
  end if
end if
!***********************************************************************
! Compute the MCSCF generalized Fock matrix and Brillouin matrix elements
!***********************************************************************
call mma_allocate(FCK,NTOT4,Label='FCK')
call mma_allocate(BM,NSXS,Label='BM')
call mma_allocate(QMat,NQ,Label='QMat') ! q-matrix(1symmblock)
call FOCK(FCK,BM,FI,FA,D,STRP,QMat,PUVX,IFINAL,CMO)
! Now FA = FI + FA. Original FA has been overwritten in FOCK routine.
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) 'FI+FA in MO-basis in sxctl (overwritten on FA)'
  write(u6,*) ' --------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if
call mma_deallocate(FCK)
call mma_deallocate(QMat)
call mma_deallocate(STRP,safe='*')
call mma_deallocate(PUVX,safe='*')

! PAM 2008: Orbital files should be updated each iteration
! for easy access in case of catastrophic failure.
if (IFINAL /= 1) then
  iShift = 0
  do ISYM=1,NSYM
    IndT = 0
    IndType(1+iShift) = NFRO(ISYM)
    IndT = IndT+NFRO(ISYM)
    IndType(2+iShift) = NISH(ISYM)
    IndT = IndT+NISH(ISYM)
    IndType(3+iShift) = NRS1(ISYM)
    IndT = IndT+NRS1(ISYM)
    IndType(4+iShift) = NRS2(ISYM)
    IndT = IndT+NRS2(ISYM)
    IndType(5+iShift) = NRS3(ISYM)
    IndT = IndT+NRS3(ISYM)
    IndType(7+iShift) = NDEL(ISYM)
    IndT = IndT+NDEL(ISYM)
    IndType(6+iShift) = NBAS(ISYM)-IndT
    iShift = iShift+7
  end do
  ! Note: This is not the final orbitals, and the orbital energies and
  ! active occupation numbers may be meaningless.
  ! There is an array with occupation numbers, so use it, even if
  ! possibly irrelevant. But put zeroes as orbital energies:
  call mma_allocate(EDUM,NTOT,Label='EDUM')
  EDUM(:) = Zero

  write(VecTyp,'(A)')
  VecTyp = '* RASSCF average (pseudo-natural) orbitals (Not final)'
  LuvvVec = 50
  LuvvVec = isfreeunit(LuvvVec)
  call WrVec('RASORB',LuvvVec,'COE',NSYM,NBAS,NBAS,CMO,OCC,EDUM,INDTYPE,VECTYP)
  call WrVec('RASORB',LuvvVec,'AI',NSYM,NBAS,NBAS,CMO,OCC,EDUM,INDTYPE,VECTYP)
  call mma_deallocate(EDUM)

# ifdef _HDF5_
  call mh5_put_dset(wfn_mocoef,CMO)
  call mh5_put_dset(wfn_occnum,OCC)
# endif
end if

if (IFINAL == 1) then
  ! If ifinal=1 , this is the last call to SXCTL for calculating
  ! the MCSCF Fock matrix for occupied orbitals, and new orbitals).
  ! First generate final orbitals:
  ! Diagonalize inactive and secondary part of FP = FI + FA.
  !         diagonal blocks of the active density matrix
  ! for RAS1, RAS2, and RAS3. All this is done in NEWORB.
  ! Finally generate orbitals and Fock matrices for CAS-PT2 runs.
  ! in FCKPT2

  ! Memory allocation and call to NEWORB and FCKPT2
  ! CMON: New molecular orbitals (NTOT2)
  ! FTR:  Temporary area for part of the Fock matrix FP (NTOT1)
  ! VEC:  EIGENVECTORS OF FTR (NO2M)
  ! SQ and WO: scratch areas

  call mma_allocate(CMON,NTOT2,Label='CMON')
  call mma_allocate(FTR,NTOT1,Label='FTR')
  call mma_allocate(VEC,NTOT2,Label='VEC')
  call mma_allocate(WO,NTOT2,Label='WO')
  call mma_allocate(SQ,NTOT2,Label='SQ')
  call mma_allocate(CMOX,NTOT2,Label='CMOX')
  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' CMO in SXCTL for IFINAL=1'
    write(u6,*) ' ---------------------'
    write(u6,*)
    ioff = 0
    do iSym=1,nSym
      iBas = nBas(iSym)
      if (iBas /= 0) then
        write(u6,*) 'Sym =',iSym
        do i=1,iBas
          write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=1,iBas)
        end do
        iOff = iOff+(iBas*iBas)
      end if
    end do
  end if

  if (iCIonly == 1) then
    IDISK = IADR15(2)
    call DDAFILE(JOBIPH,1,CMO,NTOT2,IDISK)
    call DDAFILE(JOBIPH,1,OCC,NTOT,IDISK)
#   ifdef _HDF5_
    call mh5_put_dset(wfn_mocoef,CMO)
    call mh5_put_dset(wfn_occnum,OCC)
#   endif
  else
    ! this part (TRACI) need to be changed to "TRAMPS", not yet ! Yingjin
    call NEWORB_RASSCF(CMO,CMON,FA,FTR,VEC,WO,SQ,CMOX,D,OCC)
    ! compute orbital overlap matrix
    !if (NACTEL > 0) then
    ! NN.14 Skip this when DMRG-CASSCF due to CI-vector dependency
    !if (.not. (DoDMRG .or. doBlockDMRG) .and. (NACTEL > 0)) then
    if (NACTEL > 0) then
      call mma_allocate(SMAT,NAC*NAC,Label='SMAT')
      IWAY = 1
      call OVLP(IWAY,CMO,CMON,SMAT)

      if (dodmrg) then
#       ifdef _DMRG_
#       ifdef BLUBB
        call mpsrot(mat,nac,nrs2,nsym)
#       endif
#       endif
      else if (doBlockDMRG .or. DoNECI) then
      else !CI
        iDisk = IADR15(4)
        call LUCIA_UTIL('TRACI',iDisk=iDisk,Lu=JOBIPH,Array=SMAT(:))
      end if
      call mma_deallocate(SMAT)
    else
      CIDUMMY = One
      IDISK = IADR15(4)
      call DDAFILE(JOBIPH,1,CIDUMMY,1,IDISK)
    end if
  end if

  ! IPT2 = 1 for OUTO, CANOnical option...
  if (IPT2 /= 0) call FCKPT2(CMO,CMON,FI,FA,FTR,VEC,WO,SQ,CMOX)

  call mma_deallocate(CMOX)
  call mma_deallocate(SQ)
  call mma_deallocate(WO)
  call mma_deallocate(VEC)
  call mma_deallocate(FTR)
  call mma_deallocate(CMON)

  call TIMING(CPTS,CPES,TIOS,TIOES)

  goto 9990
end if

! Memory allocation and calling sequence for SXHAM
! SXN: Normalization constants for super-CI vector
! F1 and F2: parts of the Fock matrix FP
! DIA: Occupied part of the density matrix (squared)
! SXG: The G matrix(used in sigvec)
! SXH: The H matrix( "    "   "   )
! SXHD: The diagonal of the super-CI Hamiltonian
! LDF: The matrix D*FP
! LDDIA: Diagonal of the density matrix (all elements one symmetry)

call mma_allocate(SXN,NSXS,Label='SXN')
call mma_allocate(F1,NIAIA,Label='F1')
call mma_allocate(F2,NAEAE,Label='F2')
call mma_allocate(DIA,NIAIA,Label='DIA')
call mma_allocate(SXG,NIAIA,Label='SXG')
call mma_allocate(SXH,NAOAE,Label='SXH')
call mma_allocate(SXHD,NDIMSX,Label='SXHD')
call mma_allocate(SXDF,NQ,Label='SXDF')
call mma_allocate(SXDD,MNO,Label='SXDD')

!call TRIPRT(' Dmat in MO in SXCTL bf call to SXHAM ',' ',D,NAC)
!call TRIPRT(' Pmat in MO in SXCTL bf call to SXHAM ',' ',P,NACPAR)
!call TRIPRT(' PAmat in MO in SXCTL bf call to SXHAM',' ',PA,NACPAR)
call SXHAM(D,P,PA,FA,SXN,F1,F2,DIA,SXG,SXH,SXHD,SXDF,SXDD)

call mma_deallocate(SXDD)
call mma_deallocate(SXDF)

! PAM01 Removal of certain rotations from the BLB elements.
! Some additional rotations (besides those listed in IZROT) may
! need to be suppressed. These are rotations that are (very nearly)
! redundant -- they hardly affect the wave function at all.
! All suppressed rotations can be identified because the corresponding
! diagonal elements have been set to a huge number in SXHAM.
! Use this criterion to set some BLB elements exactly =0:
do I=1,NSXS
  if (SXHD(NROOT+I) > 1.0e20_wp) BM(I) = Zero
end do

! MEMORY ALLOCATION AND CALLING SEQUENCE FOR SX DIAGONALIZATION

! CSX: The super-CI vectors
! SIGMA: The sigma vectors
! HH:  The Davidson H matrix
! CC:   "     "     egenvectors
! ENER: "     "     energies
! SC:   Scratch area
! QMat:    Davidson update vectors
! QQ:   Norm of update vectors
! OVL:  Overlap matrix

NCR = NDIMSX*NROOT*ITMAX
KMAX = ITMAX*NROOT
NCR1 = NDIMSX*NROOT*(ITMAX+1)
call mma_allocate(CSX,NCR1,Label='CSX')
call mma_allocate(SIGMA,NCR,Label='SIGMA')
NLHH = KMAX**2+KMAX
NLCC = KMAX**2
NLQ = NDIMSX*(NROOT+1)
NLOVL = ITMAX*NROOT**2
call mma_allocate(HH,NLHH,Label='HH')
call mma_allocate(CC,NLCC,Label='CC')
call mma_allocate(ENER_X,KMAX,Label='ENER_X')
call mma_allocate(SC,NDIMSX,Label='SC')
call mma_allocate(QMat,NLQ,Label='QMat')
call mma_allocate(QQ,NROOT,Label='QQ')
call mma_allocate(OVL,NLOVL,Label='OVL')

call DAVCRE(CSX,SIGMA,HH,CC,ENER_X,SXHD,SC,QMat,QQ,OVL,SXSEL,NROOT,ITMAX,NDIMSX,ITERSX,NSXS)

ESX = ENER_X(1)
call mma_deallocate(SIGMA)
call mma_deallocate(HH)
call mma_deallocate(CC)
call mma_deallocate(ENER_X)
call mma_deallocate(SC)
call mma_deallocate(QMat)
call mma_deallocate(QQ)
call mma_deallocate(OVL)
call mma_deallocate(F1)
call mma_deallocate(F2)
call mma_deallocate(SXG)
call mma_deallocate(SXH)
call mma_deallocate(SXHD)

! Renormalize the SX-coefficients

IREF = 1
LCSXI = 1+NDIMSX*(IREF-1)
IC = NROOT+LCSXI-1
XSXMAX = Zero
do I=1,NSXS
  CSX(IC+I) = SXN(I)*CSX(IC+I)/CSX(LCSXI)
  XSXMAX = max(XSXMAX,abs(CSX(IC+I)))
end do
if (IPRLEV >= DEBUG) then
  write(u6,*) 'SXCTL after DAVCRE, Renormalized SX coeffs:'
  write(u6,'(1X,8F14.6)') (CSX(IC+I),I=1,NSXS)
end if

! Step size control has been built into qune now.
!C Step length control, just for safety.
!do I=1,NSXS
!  VAL = CSX(IC+I)
!  CSX(IC+I) = VAL/(One+1.7_wp*abs(VAL))
!end do

! Intercept XSX and BM, to use (perhaps) Quasi-Newton or Line Search

!if (ITER == 1) NCALL = 0
if (ITER <= 4) NCALL = 0
if (KeyHEUR .and. (ITER > 10) .and. (mod(ITER,10) < 4)) NCALL = 0
if (doDMRG .and. (ITER <= 2)) NCALL = 0  ! YM: change 4 -> 2, for saving time
if (XSXMAX > Half) NCALL = 0
if ((NQUNE /= 0) .and. (XSXMAX < Half)) then
  call mma_allocate(VT,NSXS,Label='VT')
  call mma_allocate(VL,NSXS,Label='VL')
  call mma_allocate(XQN,NSXS,Label='XQN')
  call mma_allocate(SCR,NSXS,Label='SCR')
  call mma_allocate(V1,NSXS,Label='V1')
  call mma_allocate(V2,NSXS,Label='V2')
  CASDFT_En = Zero
  if ((KSDFT /= 'SCF') .and. (KSDFT(1:3) /= 'PAM')) call Get_dScalar('CASDFT energy',CASDFT_En)
  CASDFT_En = ECAS+CASDFT_En
  call QUNE(NCALL,CASDFT_En,BM,CSX(NROOT+LCSXI),VL,VT,XQN,SCR,V1,V2,NSXS,LUQUNE,TMIN,QNSTEP,QNUPDT,KSDFT)

  call mma_deallocate(VT)
  call mma_deallocate(VL)
  call mma_deallocate(XQN)
  call mma_deallocate(SCR)
  call mma_deallocate(V1)
  call mma_deallocate(V2)
end if

! Rotation of orbitals with exp(x) where x is obtained from
! the super-CI coefficients, with a Quasi Newton update (NQUNE=1)

! CMO:  before - old MO's           after - new MO's
! CMON: intermediate storage for new MO's (moved to CMO in ORTHO)
! X2:  work area, also in ORTHO (AO overlap matrix)
! Scr: WORK AREA

call mma_allocate(CMON,NTOT2,Label='CMON')
call mma_allocate(XMAT,NO2M,Label='XMAT')
call mma_allocate(X2,NTOT1,Label='X2')
call mma_allocate(Scr,NO2M,Label='SCR')

call ROTORB(CMO,CMON,CSX(LCSXI),XMAT,X2,SCR,THMAX,FA)

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) 'FI+FA in SXCTL after Unitary transform in ROTORB'
  write(u6,*) ' --------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if
call mma_deallocate(CMON)
call mma_deallocate(XMAT)
call mma_deallocate(X2)
call mma_deallocate(SCR)
call mma_deallocate(SXN)
call mma_deallocate(DIA)
call mma_deallocate(CSX)

IDISK = IADR15(2)
call DDAFILE(JOBIPH,1,CMO,NTOT2,IDISK)
call DDAFILE(JOBIPH,1,OCC,NTOT,IDISK)
#ifdef _HDF5_
call mh5_put_dset(wfn_mocoef,CMO)
call mh5_put_dset(wfn_occnum,OCC)
#endif
call TIMING(CPTS,CPES,TIOS,TIOES)

9990 continue
call mma_deallocate(BM)

end subroutine SXCtl
