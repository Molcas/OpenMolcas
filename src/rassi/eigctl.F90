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

!ifdef _DEBUGPRINT_
subroutine EIGCTL(PROP,OVLP,DYSAMPS,HAM,EIGVEC,ENERGY)

use Symmetry_Info, only: MUL, nIrrep
use RASSI_aux, only: iDisk_TDM, IPGLOB
use kVectors, only: e_Vector, k_Vector, nk_Vector
use rassi_global_arrays, only: JBNUM
use do_grid, only: Do_Lebedev
use nq_Grid, only: Pax
use frenkel_global_vars, only: iTyp
use rassi_data, only: NBASF, NBST, NTDMZZ
use Cntrl, only: DIPR, Do_Pol, Do_SK, Do_TMom, DoCD, DySO, EMin, iComp, IfJ2, IfJz, IPUSED, IRREP, L_Eff, LoopDivide, lSym1, &
                 lSym2, LuTDM, MLTPLT, NACTE, NPROP, nQUad, NSTATE, OSThr_DiPr, OSThr_QIPR, PNAME, PrDipVec, PrRaw, PrWeight, &
                 PTYPE, QIAll, QIPR, ReduceLoop, RSPR, RSThr, TDipMin, TMGR_Thrs, Tolerance
#ifdef _HDF5_
use Dens2HDF5, only: UpdateIdx
use mh5, only: mh5_put_dset
use Cntrl, only: IfTDM, IfTrD1
use RASSIWfn, only: wfn_SFS_Coef, wfn_SFS_Energy, WFN_SFS_TM
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Six, Nine, Ten, Half, Quart, OneHalf, Pi, auTocm, auToeV, auTofs, c_in_au, Debye
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: PROP(NSTATE,NSTATE,NPROP), OVLP(NSTATE,NSTATE), DYSAMPS(NSTATE,NSTATE), HAM(NSTATE,NSTATE), &
                 EIGVEC(NSTATE,NSTATE), ENERGY(NSTATE)
integer(kind=iwp) :: I, I2Tot, I_, I_Have_DL, I_Have_DV, I_Print_Header, iAMx, iAMxyz, iAMy, iAMz, iCar, iDiag, iDisk, idx, &
                     iEmpty, iEnd, iEnd_, IfAnyD, IfAnyM, IfAnyQ, iGo, iGrp, ii, ij, ij_, ijSO, IOFF(8), iOpt, iPrDX, iPrDxD, &
                     iPrDxM, iPrDxx, iPrDxxx, iPrDxxy, iPrDxxz, iPrDxy, iPrDxz, iPRDY, iPrDyD, iPrDyM, iPrDyx, iPrDyy, iPrDyyx, &
                     iPrDyyy, iPrDyyz, iPrDyz, iPRDZ, iPrDzD, iPrDzM, iPrDzx, iPrDzy, iPrDzz, iPrDzzx, iPrDzzy, iPrDzzz, iPrint, &
                     iProp, iPrp, iPrQxx, iPrQxy, iPrQxz, iPrQyy, iPrQyz, iPrQzz, IPRTMOM(12), iQuad, iSet, iStart_, iState, &
                     iSy12, iSy34, iTol, iType, iVec, iVec_, J, J_, jEnd_, jGrp, jj, Job1, Job2, Job3, Job4, jSet, jSTart, &
                     jStart_, jState, k, K_, kEnd, kSTA, kState, L, L_, LNCNT, lOsc_Strength, lPos, lState, lSym3, lSym4, LuT1, &
                     Mask, Mask34, MaxGrp1, MaxGrp2, MPLET1, MPLET2, mState, nActe1, nActe2, nDiff, nGroup1, nGroup2, nHH, NIP, &
                     nLST, nMax2, nScr, nSets, nTmp, nVec, SECORD(4)
real(kind=wp) :: A, AFactor, ANG, Ax, Ay, Az, COMPARE, DMax, DSZ, Dx, Dx2, Dxx, Dxx2, DxxDyy, DxxDzz, DxxxDx, DxxyDy, DxxzDz, Dxy, &
                 Dxy2, DxyDz, Dxz, Dxz2, DxzDy, Dy, Dy2, DYSAMPS2(NSTATE,NSTATE), DysThr, DyxDz, Dyy, Dyy2, DyyDzz, DyyxDx, &
                 DyyyDy, DyyzDz, Dyz, Dyz2, DyzDx, Dz, Dz2, DzxDy, DzyDx, Dzz, Dzz2, DzzxDx, DzzyDy, DzzzDz, E0, E1, E2, E3, &
                 EDiff, EDiff2, EDiff3, EDiff_, EffL, EffM, EI, epsH, epsS, eRMS, EV, EVLim, EVMax, Ex, F, F_Check, F_Temp, FMax, &
                 Fx, Fxx, FxxFyy, FxxFzz, Fxxx, Fxxy, Fxxz, Fxy, Fxz, Fy, Fyx, Fyy, FyyFzz, Fyyx, Fyyy, Fyyz, Fyz, Fz, Fzx, FZY, &
                 Fzz, Fzzx, Fzzy, Fzzz, OSthr, OSthr2, R, R_Check, R_Temp, RefEne, RKNorm, RNG, RNorm, Rtensor(6), Rxx, Rxxy, &
                 Rxxz, Rxy, Rxyx, Rxyy, Rxyz, Rxz, Rxzx, Rxzy, Rxzz, Ryx, Ryy, Ryyx, Ryyz, Ryz, Ryzx, Ryzy, RYZZ, Rzx, Rzy, Rzz, &
                 Rzzx, Rzzy, Tau, TCPU1, TCPU2, ThrS, TM1, TM2, TM3, TM_2, TM_C(3), TM_I(3), TM_R(3), Tmp, TWall1, TWall2, UK(3), &
                 V2Sum, Wavevector(3), Weight, X
logical(kind=iwp) :: Diagonal, TMOgroup
character(len=100) :: line
character(len=60) :: FMTLINE
character(len=13) :: filnam
character(len=8) :: LABEL
integer(kind=iwp), allocatable :: ILST(:), IndexE(:), LIST(:), STACK(:), TMOgrp1(:), TMOgrp2(:)
real(kind=wp), allocatable :: Aux(:,:), DL(:), DV(:), ESFS(:), HH(:), HSQ(:), IP(:), L2(:), L2DIA(:), M2DIA(:), OscStr(:,:), &
                              pol_Vector(:,:), RAW(:,:,:), Rquad(:,:), SCR(:,:), SCR1(:), SS(:), TDMZZ(:), TOT2K(:,:), TSDMZZ(:), &
                              UU(:), VLST(:), WDMZZ(:)
#ifdef _HDF5_
integer(kind=iwp) :: ijSF, ip_kVector, ip_TMI, ip_TMR, ip_W, nData, nIJ
real(kind=wp), allocatable :: Storage(:,:,:,:)
#endif
real(kind=wp), parameter :: AU2REDR = 200.0_wp*Debye, DLT = 1.0e-18_wp, ONEOVER10C = One/(Ten*c_in_au**2), &
                            ONEOVER30C = ONEOVER10C/Three, ONEOVER6C2 = One/(Six*c_in_au**2), ONEOVER9C2 = One/(Nine*c_in_au**2), &
                            Two3rds = Two/Three, TWOOVERM45C = -Two/(45.0_wp*c_in_au**2)
integer(kind=iwp), external :: cho_x_gettol, IsFreeUnit
real(kind=wp), external :: DDOt_

! Bruno, DYSAMPS2 is used for printing out the pure norm
! of the Dyson vectors.
! DYSAMPS remains basis of the SF eigen-states to the basis
! of the original SF states.
DYSAMPS2 = DYSAMPS

DIAGONAL = .true.

#ifdef _DEBUGPRINT_
write(u6,*) 'BLUBB start of eigctl: debug print of property matrix'
do istate=1,nstate
  do jstate=1,nstate
    do IPROP=1,NPROP
      if (abs(prop(istate,jstate,iprop)) > 1.0e-14_wp) &
        write(u6,*) 'prop(',istate,',',jstate,',',iprop,') = ',prop(istate,jstate,iprop)
    end do
  end do
end do
#endif

! DIAGONALIZE SCALAR HAMILTONIAN.

! Initialize eigenvector array.
EIGVEC(:,:) = Zero
! NOTE: It is imperative that we do not mix, or change order, of
! states with different nr of electrons, spin, or symmetry. It is
! assumed in some subsequent parts of the program that these
! remain good properties of the wave functions, with values as
! listed in various tables in common /CNTRL/.
! So it is worth the extra inconvenience to construct an outer
! loop over sets of interacting wave functions.
! Make a list of interacting sets of states:
call mma_allocate(LIST,NSTATE,Label='LIST')
LIST(:) = 0
ISET = 0
do I=1,NSTATE
  if (LIST(I) > 0) cycle
  ISET = ISET+1
  LIST(I) = ISET
  JOB1 = JBNUM(I)
  NACTE1 = NACTE(JOB1)
  MPLET1 = MLTPLT(JOB1)
  LSYM1 = IRREP(JOB1)
  do J=I+1,NSTATE
    if (LIST(J) > 0) cycle
    JOB2 = JBNUM(J)
    NACTE2 = NACTE(JOB2)
    if (NACTE2 /= NACTE1) cycle
    MPLET2 = MLTPLT(JOB2)
    if (MPLET2 /= MPLET1) cycle
    LSYM2 = IRREP(JOB2)
    if (LSYM2 /= LSYM1) cycle
    LIST(J) = ISET
  end do
end do
NSETS = ISET
!TEST write(u6,*)' EIGCTL. There are NSETS sets of interacting states.'
!TEST write(u6,'(1x,a,i3)')' where NSETS=',NSETS
!TEST write(u6,*)' The LIST array:'
!TEST write(u6,'(1x,20i3)') (LIST(I),I=1,NSTATE)

NHH = (NSTATE*(NSTATE+1))/2
call mma_allocate(HH,NHH,Label='HH')
call mma_allocate(HSQ,NSTATE**2,Label='HSQ')
call mma_allocate(SS,NHH,Label='SS')
call mma_allocate(UU,NSTATE**2,Label='UU')
call mma_allocate(SCR1,NSTATE**2,Label='SCR1')
SCR1(:) = Zero
call mma_allocate(STACK,NSTATE,Label='STACK')
! Loop over the sets:
do ISET=1,NSETS
  ! Stack up the states belonging to this set:
  MSTATE = 0
  do I=1,NSTATE
    JSET = LIST(I)
    if (JSET == ISET) then
      MSTATE = MSTATE+1
      STACK(MSTATE) = I
    end if
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'BLUBB DEBUG print of Hamiltonian and overlap'
# endif

  ! 1. PUT UNIT MATRIX INTO UU
  UU(:) = Zero
  call DCOPY_(MSTATE,[One],0,UU,MSTATE+1)
  ! 2. COPY OVERLAP MATRIX INTO TRIANGULAR STORAGE,
  !    and Hamiltonian into square storage:
  call DCOPY_(NSTATE**2,[Zero],0,HSQ,1)
  IJ = 0
  do II=1,MSTATE
    I = STACK(II)
    do JJ=1,II
      J = STACK(JJ)
      IJ = IJ+1
      if ((I /= J) .and. ((abs(ovlp(i,j)) > 1.0e-9_wp) .or. (abs(ham(i,j)) > 1.0e-9_wp))) DIAGONAL = .false.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'overlap     for i,j',i,j,ovlp(i,j)
      write(u6,*) 'Hamiltonian for i,j',i,j,HAM(i,j)
#     endif
      SS(IJ) = OVLP(I,J)
      HSQ(II+MSTATE*(JJ-1)) = HAM(I,J)
      HSQ(JJ+MSTATE*(II-1)) = HAM(I,J)
    end do
  end do
  ! 3. SPECTRAL DECOMPOSITION OF OVERLAP MATRIX:
  call Jacob(SS,UU,MSTATE,MSTATE)
  II = 0
  do I=1,MSTATE
    II = II+I
    X = One/sqrt(max(5.0e-15_wp,SS(II)))
    do K=1,MSTATE
      LPOS = K+MSTATE*(I-1)
      UU(LPOS) = X*UU(LPOS)
    end do
  end do

  if (.not. diagonal) then
    ! 4. TRANSFORM HAMILTON MATRIX.
    call DGEMM_('N','N',MSTATE,MSTATE,MSTATE,One,HSQ,MSTATE,UU,MSTATE,Zero,SCR1,MSTATE)
    call DGEMM_('T','N',MSTATE,MSTATE,MSTATE,One,UU,MSTATE,SCR1,MSTATE,Zero,HSQ,MSTATE)

    ! 5. DIAGONALIZE HAMILTONIAN.
    IJ = 0
    do I=1,MSTATE
      do J=1,I
        IJ = IJ+1
        HH(IJ) = HSQ(I+MSTATE*(J-1))
      end do
    end do

    call Jacob(HH,UU,MSTATE,MSTATE)
    call SortDiag(HH,UU,MSTATE,MSTATE)

    IDIAG = 0
    do II=1,MSTATE
      IDIAG = IDIAG+II
      I = STACK(II)
      ENERGY(I) = HH(IDIAG)
      do JJ=1,MSTATE
        J = STACK(JJ)
        EIGVEC(I,J) = UU(II+MSTATE*(JJ-1))
      end do
    end do
  else
    !if diagonal
    do II=1,MSTATE
      I = STACK(II)
      ENERGY(I) = HAM(I,I)
      do JJ=1,MSTATE
        J = STACK(JJ)
        EIGVEC(I,J) = Zero
      end do
      EIGVEC(I,I) = One
    end do
  end if

  !UNGUR
  !   Correct for diagonal energies in case of orbital degeneracy:
  !   Convention: two energies are considered degenerate if their energy difference is
  !               lower than 1.0e-4 cm-1
  TMP = Zero
  IDIAG = 0
  do II=1,MSTATE
    I = STACK(II)
    TMP = ENERGY(I)
    do JJ=1,MSTATE
      J = STACK(JJ)
      if (I == J) cycle
      if (abs(ENERGY(J)-TMP)*auTocm < 1.0e-4_wp) ENERGY(J) = TMP
    end do
  end do

  IDIAG = 0
  do II=1,MSTATE
    IDIAG = IDIAG+II
    I = STACK(II)
    HH(IDIAG) = ENERGY(I)
  end do
  ! End of loop over sets.
end do
! Morgane Vacher 02/17 - Fix the "arbitrary" sign of
! the eigenvectors such that the largest coefficient
! is positive. This is to avoid spurious changes of
! sign of the SFS with respect to the original ones,
! especially for already diagonal Hamiltonian matrix.
do i=1,nstate
  j = maxloc(abs(eigvec(:,i)),1)
  if (eigvec(i,j) < Zero) eigvec(:,i) = -eigvec(:,i)
end do
call mma_deallocate(HH)
call mma_deallocate(SS)
call mma_deallocate(UU)
call mma_deallocate(HSQ)
call mma_deallocate(STACK)
call mma_deallocate(LIST)

#ifdef _HDF5_
call mh5_put_dset(wfn_sfs_energy,ENERGY)
call mh5_put_dset(wfn_sfs_coef,EIGVEC)
#endif

if (IPGLOB >= 1) then
  write(filnam,'(A,I1)') 'stE',iTyp
  LuT1 = isFreeUnit(11)
  call molcas_open(LuT1,filnam)
  do istate=1,nstate
    write(LuT1,*) energy(istate)
  end do
  close(LuT1)
  do istate=1,nstate
    call PrintResult(u6,'(6x,A,I5,5X,A,F23.14)','RASSI State',ISTATE,'Total energy:',ENERGY(ISTATE),1)
  end do
end if

! Put energies onto info file for automatic verification runs:
!PAM06 Added error estimate, based on independent errors for all
! components of H and S in original RASSCF wave function basis:
EPSS = 5.0e-11_wp
EPSH = max(5.0e-10_wp,abs(ENERGY(1))*EPSS)
IDX = 100
do I=1,NSTATE
  EI = ENERGY(I)*EPSS
  V2SUM = sum(EIGVEC(:,I)**2)
  ERMS = sqrt(EPSH**2+EI**2)*V2SUM
  IDX = min(IDX,int(-log10(ERMS)))
end do
iTol = cho_x_gettol(IDX) ! reset thr iff Cholesky
call Add_Info('E_RASSI',ENERGY,NSTATE,iTol)

! To handle extreme cases of large energies/small energy differences
! all TOTAL energies will undergo a universal constant shift:
EMIN = minval(ENERGY(:))
ENERGY(:) = ENERGY(:)-EMIN

! Experimental addition: Effective L and/or M quantum numbers.

! Identify which properties are angular moment matrix elements:
IAMX = 0
IAMY = 0
IAMZ = 0
do IPROP=1,NPROP
  if (PNAME(IPROP)(1:6) == 'ANGMOM') then
    if (ICOMP(IPROP) == 1) IAMX = IPROP
    if (ICOMP(IPROP) == 2) IAMY = IPROP
    if (ICOMP(IPROP) == 3) IAMZ = IPROP
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! The matrix elements are actually for i*Lx, etc.
! Now form matrix elements of L**2 assuming closure:
! i.e. assume L**2 = -(iLx)*(iLx)-(iLy)*(iLy)-(iLz)*(iLz)
! within the basis formed by the states.
IAMXYZ = 0
if (IAMZ > 0) then
  call mma_allocate(L2,NSTATE**2,Label='L2')
  call mma_allocate(M2DIA,NSTATE,Label='M2DIA')
  call DGEMM_('N','N',NSTATE,NSTATE,NSTATE,-One,PROP(1,1,IAMZ),NSTATE,PROP(1,1,IAMZ),NSTATE,Zero,L2,NSTATE)
  call DCOPY_(NSTATE,L2,(NSTATE+1),M2DIA,1)
  if ((IAMX > 0) .and. (IAMY > 0)) then
    IAMXYZ = 1
    call DGEMM_('N','N',NSTATE,NSTATE,NSTATE,-One,PROP(1,1,IAMX),NSTATE,PROP(1,1,IAMX),NSTATE,One,L2,NSTATE)
    call DGEMM_('N','N',NSTATE,NSTATE,NSTATE,-One,PROP(1,1,IAMY),NSTATE,PROP(1,1,IAMY),NSTATE,One,L2,NSTATE)
    call mma_allocate(L2DIA,NSTATE,Label='L2DIA')
    call DCOPY_(NSTATE,L2,(NSTATE+1),L2DIA,1)
  end if
  call mma_deallocate(L2)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Sort the states energywise

call mma_Allocate(IndexE,nState,Label='IndexE')
IndexE(:) = [(iState,iState=1,nState)]
do iState=1,nState-1
  EX = ENERGY(IndexE(iState))

  kState = iState
  do jState=iState+1,nState
    if (ENERGY(IndexE(jState)) < EX) then
      kState = jState
      EX = ENERGY(IndexE(jState))
    end if
  end do
  if (kState /= iState) then
    lState = IndexE(iState)
    IndexE(iState) = IndexE(kState)
    IndexE(kState) = lState
  end if
end do
#ifdef _HDF5_
if (IFTRD1 .or. IFTDM) call UpdateIdx(IndexE,0)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! REPORT ON SECULAR EQUATION RESULT:
call MMA_ALLOCATE(ESFS,NSTATE)
if (IPGLOB >= 1) then
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A,34X,A,34X,A)') '*','       Spin-free section      ','*'
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A,17X,A,17X,A)') '*','Note: index according to input order, order according to energy.','*'
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,*) ' SPIN-FREE ENERGIES:'
  write(u6,'(1X,A,F22.10,A1)') ' (Shifted by EMIN (a.u.) =',EMIN,')'
  write(u6,*)
  if ((IFJ2 /= 0) .and. (IAMXYZ > 0)) then
    if ((IFJZ /= 0) .and. (IAMZ > 0)) then
      write(u6,*) 'SF State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)      L_eff   Abs_M'
    else
      write(u6,*) 'SF State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)      L_eff'
    end if
  else
    if ((IFJZ /= 0) .and. (IAMZ > 0)) then
      write(u6,*) 'SF State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)      Abs_M'
    else
      write(u6,*) 'SF State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)'
    end if
  end if
  write(u6,*)

  E0 = ENERGY(IndexE(1))
  do kSTATE=1,NSTATE
    iState = IndexE(kState)
    E1 = ENERGY(ISTATE)
    E2 = auToeV*(E1-E0)
    E3 = auTocm*(E1-E0)

    if ((IFJ2 /= 0) .and. (IAMXYZ > 0)) then
      if ((IFJZ /= 0) .and. (IAMZ > 0)) then
        EFFL = sqrt(max(5.0e-13_wp,Quart+L2DIA(ISTATE)))-Half
        EFFM = sqrt(max(5.0e-13_wp,M2DIA(ISTATE)))
        FMTLINE = '(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4,6X,F6.1,2X,F6.1)'
        write(u6,FMTLINE) ISTATE,E1,E2,E3,EFFL,EFFM
      else
        EFFL = sqrt(max(5.0e-13_wp,Quart+L2DIA(ISTATE)))-Half
        FMTLINE = '(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4,6X,F6.1)'
        write(u6,FMTLINE) ISTATE,E1,E2,E3,EFFL
      end if
    else
      if ((IFJZ /= 0) .and. (IAMZ > 0)) then
        EFFM = sqrt(max(5.0e-13_wp,M2DIA(ISTATE)))
        FMTLINE = '(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4,6X,F6.1)'
        write(u6,FMTLINE) ISTATE,E1,E2,E3,EFFM
      else
        FMTLINE = '(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4)'
        write(u6,FMTLINE) ISTATE,E1,E2,E3
      end if
    end if
    ESFS(ISTATE) = E3

  end do

end if
if (IAMZ > 0) call mma_deallocate(M2DIA)
if (IAMXYZ > 0) call mma_deallocate(L2DIA)
! LU: save esfs array
call Put_dArray('ESFS_SINGLE',ESFS,NSTATE)
call Put_dArray('ESFS_SINGLEAU',(ENERGY+EMIN),NSTATE)
call MMA_DEALLOCATE(ESFS)

if ((IPGLOB >= 3) .or. (.not. diagonal)) then
  write(u6,*)
  write(u6,*) '  Spin-free eigenstates in basis of input states:'
  write(u6,*) '  -----------------------------------------------'
  write(u6,*)
  if (IPGLOB >= 3) then
    do L=1,NSTATE
      I = IndexE(L)
      write(u6,'(5X,A,I5,A,F18.10)') 'Eigenstate No.',I,' energy=',ENERGY(I)+EMIN
      write(u6,'(5X,5F15.7)') (EIGVEC(K,I),K=1,NSTATE)
    end do
  end if
  call mma_allocate(ILST,NSTATE,Label='ILST')
  call mma_allocate(VLST,NSTATE,Label='VLST')
  do L=1,NSTATE
    I = IndexE(L)
    write(u6,'(5X,A,I5,A,F18.10)') 'Eigenstate No.',I,' energy=',ENERGY(I)+EMIN
    EVMAX = maxval(abs(EIGVEC(:,I)))
    EVLIM = 0.1_wp*EVMAX
    NLST = 0
    do K=1,NSTATE
      EV = EIGVEC(IndexE(K),I)
      if (abs(EV) >= EVLIM) then
        NLST = NLST+1
        VLST(NLST) = EV
        ILST(NLST) = IndexE(K)
      end if
    end do
    do KSTA=1,NLST,6
      KEND = min(NLST,KSTA+4)
      write(Line,'(5X,5(I5,F12.6))') (ILST(K),VLST(K),K=KSTA,KEND)
      call NORMAL(Line)
      write(u6,*) Line
    end do
    write(u6,*)
  end do
  call mma_deallocate(ILST)
  call mma_deallocate(VLST)
  if (IPGLOB >= 3) then
    write(u6,*)
    write(u6,*) ' THE INPUT RASSCF STATES REEXPRESSED IN EIGENSTATES:'
    write(u6,*)
    do L=1,NSTATE
      I = IndexE(L)
      call DGEMM_('T','N',NSTATE,NSTATE,NSTATE,One,EIGVEC,NSTATE,OVLP,NSTATE,Zero,SCR1,NSTATE)
      write(u6,'(A,I5)') ' INPUT STATE NR.:',I
      write(u6,*) ' OVERLAP WITH THE EIGENSTATES:'
      write(u6,'(5(1X,F15.7))') (SCR1(IndexE(K)+NSTATE*(I-1)),K=1,NSTATE)
      write(u6,*)
    end do
  end if
end if
!                                                                      !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!                                                                      !
!                                                                      !
!     TRANSFORM AND PRINT OUT PROPERTY MATRICES:                       !
!                                                                      !
!     The matrix elements of Prop refers to the original SF basis. We  !
!     now transform these to the basis of the eigenvectors of the SF   !
!     states. Note, for the "exact" operator of the transition moments !
!     we will have a computation of the TDM of the SF states on the    !
!     fly. To account for this transformation we will have to transform!
!     the coefficients of the SO states from a basis of the SF eigen-  !
!     states to the basis of the original SF states.                   !

do IPRP=1,NPROP
  call DGEMM_('N','N',NSTATE,NSTATE,NSTATE,One,PROP(1,1,IPRP),NSTATE,EIGVEC,NSTATE,Zero,SCR1,NSTATE)
  call DGEMM_('T','N',NSTATE,NSTATE,NSTATE,One,EIGVEC,NSTATE,SCR1,NSTATE,Zero,PROP(1,1,IPRP),NSTATE)
end do

! And the same for the Dyson amplitudes
if (DYSO) then
  call DGEMM_('N','N',NSTATE,NSTATE,NSTATE,One,DYSAMPS,NSTATE,EIGVEC,NSTATE,Zero,SCR1,NSTATE)
  call DGEMM_('T','N',NSTATE,NSTATE,NSTATE,One,EIGVEC,NSTATE,SCR1,NSTATE,Zero,DYSAMPS,NSTATE)
end if
!                                                                      !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!                                                                      !

call mma_deallocate(SCR1)

! Initial setup for both dipole, quadrupole etc. and exact operator

! There are debug statements thoughout - look for ZVAL
! If you want to debug in length gauge then comment out velocity dipole

!ZVAL(1) = One ! Simulation of moving the origin along Z
!ZVAL(2) = Two
!ZVAL(3) = Three
!ZVAL(4) = Five
!ZVAL(5) = Seven
!ZVAL(6) = Ten
!ZVAL(7) = 15.0_wp
!ZVAL(8) = 20.0_wp
!ZVAL(9) = 25.0_wp

OSTHR = 1.0e-5_wp
if (DIPR) OSTHR = OSTHR_DIPR
if (DIPR) write(u6,30) 'Dipole printing threshold changed to ',OSTHR
! this is to ensure that the total transition strength is non-zero
! Negative transitions strengths can occur for quadrupole transitions
! due to the truncation of the Taylor expansion.
if (QIPR) then
  OSTHR = OSTHR_QIPR
  write(u6,49) 'Printing threshold changed to ',OSTHR,' since quadrupole printing threshold is given'
end if
OSTHR2 = 1.0e-5_wp
if (QIPR) then
  OSTHR2 = OSTHR_QIPR
  write(u6,30) 'Quadrupole printing threshold changed to ',OSTHR2
end if
if (QIALL) write(u6,*) 'Will write all quadrupole contributions'

! Rotatory strength threshold
if (RSPR) then
  write(u6,30) 'Rotatory strength printing threshold changed to ',RSTHR
else
  RSTHR = 1.0e-7_wp ! Default
end if

! Reducing the loop over states - good for X-rays
! At the moment memory is not reduced

if (REDUCELOOP) then
  IEND = LOOPDIVIDE
  JSTART = LOOPDIVIDE+1
else
  IEND = NSTATE
  JSTART = 1
end if

! AFACTOR = 2*pi*e^2*E_h^2 / eps_0*m_e*c^3*h^2
! numerically: 2/c^3 (in a.u. of time ^ -1)
AFACTOR = Two/c_in_au**3/(auTofs*1.0e-15_wp)

if (IPGLOB <= 0) then
  call FinishUp()
  return
end if

! CALCULATION OF THE DIPOLE TRANSITION STRENGTHS

! Initialize arrays for indentifying problematic transitions
! These stores all dipole oscillator strengths in
! length and velocity gauge for a later comparison.

call mma_allocate(DL,NSTATE**2,Label='DL')
call mma_allocate(DV,NSTATE**2,Label='DV')
DL(:) = Zero
DV(:) = Zero
I_HAVE_DL = 0
I_HAVE_DV = 0

if (IPGLOB >= 1) then

  IPRDX = 0
  IPRDY = 0
  IPRDZ = 0
  IFANYD = 0
  do IPROP=1,NPROP
    if (PTYPE(IPROP)(5:8) /= 'SING') cycle
    if (IPUSED(IPROP) /= 0) then
      if (PNAME(IPROP) == 'MLTPL  1') then
        IFANYD = 1
        if (ICOMP(IPROP) == 1) IPRDX = IPROP
        if (ICOMP(IPROP) == 2) IPRDY = IPROP
        if (ICOMP(IPROP) == 3) IPRDZ = IPROP
      end if
    end if
  end do

  if (IFANYD /= 0) then
    write(u6,*)
    call CollapseOutput(1,'Dipole transition strengths (spin-free states):')
    write(u6,'(3X,A)') '-----------------------------------------------'
    if (OSTHR > Zero) write(u6,30) 'for osc. strength at least',OSTHR
    write(u6,*)

    ! Printout the osc. strength in 3 dimensions into a file
    ! Should be if something happen
    losc_strength = 20
    losc_strength = isFreeUnit(losc_strength)
    call Molcas_Open(losc_strength,'osc_strength.au')

    if (Do_SK) then
      nVec = nk_Vector
    else
      nVec = 1
    end if

    do iVec=1,nVec

      LNCNT = 0
      FMAX = Zero
      do K_=1,IEND
        I = IndexE(K_)
        do L_=JSTART,NSTATE
          J = IndexE(L_)
          IJ = I+NSTATE*(J-1)
          EDIFF = ENERGY(J)-ENERGY(I)

          if (EDIFF > Zero) then
            DX = Zero
            DY = Zero
            DZ = Zero
            if (IPRDX > 0) DX = PROP(J,I,IPRDX)
            if (IPRDY > 0) DY = PROP(J,I,IPRDY)
            if (IPRDZ > 0) DZ = PROP(J,I,IPRDZ)
            if (Do_SK) then
              tmp = DX*k_vector(1,iVec)+DY*k_vector(2,iVec)+DZ*k_vector(3,iVec)
              DX = DX-tmp*k_vector(1,iVec)
              DY = DY-tmp*k_vector(2,iVec)
              DZ = DZ-tmp*k_vector(3,iVec)
            end if
            DX2 = DX**2
            DY2 = DY**2
            DZ2 = DZ**2
            FX = Two3rds*EDIFF*(DX2)
            FY = Two3rds*EDIFF*(DY2)
            FZ = Two3rds*EDIFF*(DZ2)
            F = FX+FY+FZ
            FMAX = max(F,FMAX)
            AX = (AFACTOR*EDIFF**2)*FX
            AY = (AFACTOR*EDIFF**2)*FY
            AZ = (AFACTOR*EDIFF**2)*FZ
            A = (AFACTOR*EDIFF**2)*F
            if (F >= OSTHR) then
              if (LNCNT == 0) then
                if (Do_SK) then
                  write(u6,*)
                  write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3)
                  write(u6,'(4x,a)') 'The light is assumed to be unpolarized.'
                  write(u6,*)
                end if
                write(u6,31) 'From','To','Osc. strength','Einstein coefficients Ax, Ay, Az (sec-1)   ','Total A (sec-1)'
                write(u6,32)
                write(losc_strength,34) 'From','To','Osc. strength','Fx','Fy','Fz','(a.u.)'
                write(losc_strength,32)
              end if
              LNCNT = LNCNT+1
              write(u6,33) I,J,F,AX,AY,AZ,A
              write(losc_strength,33) I,J,F,Fx,Fy,Fz
            end if
            ! Store dipole oscillator strength
            DL(IJ) = F

            if (F > One) then
              k = int(log10(F))+1
              F = F/(Ten**k)
            end if
            call Add_Info('TMS(SF,Len)',[F],1,6)
          end if
        end do
      end do
      if (LNCNT == 0) then
        write(u6,*) ' ( Max oscillator strength is only ',FMAX,')'
      else
        write(u6,32)
        write(losc_strength,32)
      end if

    end do ! iVec

    close(losc_strength)
    call CollapseOutput(0,'Dipole transition strengths (spin-free states):')
    I_HAVE_DL = 1
  end if

  ! Keywords for printing transition dipole vectors
  ! PRDIPVEC TDIPMIN
  if (PRDIPVEC .and. (NSTATE > 1) .and. (IFANYD /= 0)) then
    write(filnam,'(A,I1)') 'dip_vec',iTyp
    LuT1 = isFreeUnit(11)
    call molcas_open(LuT1,filnam)

    write(u6,*)
    call CollapseOutput(1,'Dipole transition vectors (spin-free states):')
    write(u6,'(3X,A)') '---------------------------------------------'
    if (TDIPMIN > Zero) write(u6,30) 'for vector sizes at least',TDIPMIN
    write(u6,*)

    if (Do_SK) then
      nVec = nk_Vector
    else
      nVec = 1
    end if

    do iVec=1,nVec

      LNCNT = 0
      DMAX = Zero
      do K_=1,IEND
        I = IndexE(K_)
        do L_=JSTART,NSTATE
          J = IndexE(L_)
          EDIFF = ENERGY(J)-ENERGY(I)
          if (EDIFF > Zero) then
            DX = Zero
            DY = Zero
            DZ = Zero
            if (IPRDX > 0) DX = PROP(J,I,IPRDX)
            if (IPRDY > 0) DY = PROP(J,I,IPRDY)
            if (IPRDZ > 0) DZ = PROP(J,I,IPRDZ)
            if (Do_SK) then
              tmp = DX*k_vector(1,iVec)+DY*k_vector(2,iVec)+DZ*k_vector(3,iVec)
              DX = DX-tmp*k_vector(1,iVec)
              DY = DY-tmp*k_vector(2,iVec)
              DZ = DZ-tmp*k_vector(3,iVec)
            end if
            DSZ = sqrt(DX**2+DY**2+DZ**2)
            DMAX = max(DSZ,DMAX)
            if (DSZ >= TDIPMIN) then
              if (LNCNT == 0) then
                if (Do_SK) then
                  write(u6,*)
                  write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3)
                  write(u6,'(4x,a)') 'The light is assumed to be unpolarized.'
                  write(u6,*)
                end if
                write(u6,34) 'From','To','Dx','Dy','Dz','Total D (a.u.)'
                write(u6,42)
              end if
              LNCNT = LNCNT+1
              write(u6,33) I,J,DX,DY,DZ,DSZ
              write(LuT1,222) I,J,DX,DY,DZ
            end if
            call Add_Info('TVC(SF,Len)',[DSZ],1,6)
          end if
        end do
      end do

      if (LNCNT == 0) then
        write(u6,*) ' ( Max transition dipole is only ',DMAX,')'
      else
        write(u6,32)
      end if
      call CollapseOutput(0,'Dipole transition vectors (spin-free states):')

    end do ! iVec

  end if
  close(LuT1)
end if

! Transition moments computed with the velocity operator.

IPRDX = 0
IPRDY = 0
IPRDZ = 0
IFANYD = 0
do IPROP=1,NPROP
  if (IPUSED(IPROP) /= 0) then
    if (PNAME(IPROP) == 'VELOCITY') then
      IFANYD = 1
      if (ICOMP(IPROP) == 1) IPRDX = IPROP
      if (ICOMP(IPROP) == 2) IPRDY = IPROP
      if (ICOMP(IPROP) == 3) IPRDZ = IPROP
    end if
  end if
end do

if (IFANYD /= 0) then
  write(u6,*)
  call CollapseOutput(1,'Velocity transition strengths (spin-free states):')
  write(u6,'(3X,A)') '-------------------------------------------------'
  if (OSTHR > Zero) write(u6,30) 'for osc. strength at least',OSTHR
  write(u6,*)

  if (Do_SK) then
    nVec = nk_Vector
  else
    nVec = 1
  end if

  do iVec=1,nVec

    LNCNT = 0
    FMAX = Zero
    do K_=1,IEND
      I = IndexE(K_)
      do L_=JSTART,NSTATE
        J = IndexE(L_)
        IJ = I+NSTATE*(J-1)
        EDIFF = ENERGY(J)-ENERGY(I)
        if (EDIFF > Zero) then
          DX = Zero
          DY = Zero
          DZ = Zero
          if (IPRDX > 0) DX = PROP(J,I,IPRDX)
          if (IPRDY > 0) DY = PROP(J,I,IPRDY)
          if (IPRDZ > 0) DZ = PROP(J,I,IPRDZ)
          if (Do_SK) then
            tmp = DX*k_vector(1,iVec)+DY*k_vector(2,iVec)+DZ*k_vector(3,iVec)
            DX = DX-tmp*k_vector(1,iVec)
            DY = DY-tmp*k_vector(2,iVec)
            DZ = DZ-tmp*k_vector(3,iVec)
          end if
          DX2 = DX**2
          DY2 = DY**2
          DZ2 = DZ**2
          FX = Two3rds*(DX2)/EDIFF
          FY = Two3rds*(DY2)/EDIFF
          FZ = Two3rds*(DZ2)/EDIFF
          F = FX+FY+FZ
          FMAX = max(F,FMAX)
          AX = (AFACTOR*EDIFF**2)*FX
          AY = (AFACTOR*EDIFF**2)*FY
          AZ = (AFACTOR*EDIFF**2)*FZ
          A = (AFACTOR*EDIFF**2)*F
          if (F >= OSTHR) then
            if (LNCNT == 0) then
              if (Do_SK) then
                write(u6,*)
                write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,ivec),k=1,3)
                write(u6,'(4x,a)') 'The light is assumed to be unpolarized.'
                write(u6,*)
              end if
              write(u6,31) 'From','To','Osc. strength','Einstein coefficients Ax, Ay, Az (sec-1)   ','Total A (sec-1)'
              write(u6,32)
            end if
            LNCNT = 1
            write(u6,33) I,J,F,AX,AY,AZ,A
            ! Store dipole oscillator strength
            DV(IJ) = F
          end if
          call Add_Info('TMS(SF,Vel)',[F],1,6)
        end if
      end do
    end do
    if (LNCNT == 0) then
      write(u6,*) ' ( Max oscillator strength is only ',FMAX,')'
    else
      write(u6,32)
    end if

  end do ! iVec
  call CollapseOutput(0,'Velocity transition strengths (spin-free states):')
  write(u6,*)
  I_HAVE_DV = 1
end if

! Compare oscillator strengths in length and velocity gauge
! All differences in oscillator strengths above the tolerance
! of 0.1 (10 percent) will be printed.

if ((I_HAVE_DL == 1) .and. (I_HAVE_DV == 1)) then
  call CollapseOutput(1,'Length and velocity gauge comparison (spin-free states):')
  write(u6,'(3X,A)') '--------------------------------------------------------'

  ! I guess that I have to explain it when I print a warning

  write(u6,*)
  write(u6,*) '--------------------------------------------------'
  write(u6,*) 'A comparison between the dipole oscillator strengths in'
  write(u6,*) 'length and velocity gauge will be performed'
  write(u6,*)
  write(u6,49) 'All dipole oscillator differences above the tolerance of ',TOLERANCE,' will be printed'
  write(u6,*)
  write(u6,*) 'Due to basis set deficiency these oscillator may be problematic'
  write(u6,*)
  write(u6,*) 'The tolerance is defined as ABS(1-O_l/O_v)'
  write(u6,*) 'O_l : dipole oscillator strength in length gauge'
  write(u6,*) 'O_p : dipole oscillator strength in velocity gauge'
  write(u6,*) '--------------------------------------------------'

  I_PRINT_HEADER = 0
  do K_=1,IEND
    I = IndexE(K_)
    do L_=JSTART,NSTATE
      J = IndexE(L_)
      IJ = I+NSTATE*(J-1)
      EDIFF = ENERGY(J)-ENERGY(I)
      if ((JSTART == 1) .and. (EDIFF < Zero)) cycle

      COMPARE = Zero
      if ((DL(IJ) >= OSTHR+dlt) .and. (DV(IJ) >= OSTHR+dlt)) then
        COMPARE = abs(1-DL(IJ)/DV(IJ))
      else if ((DL(IJ) >= OSTHR+dlt) .and. (DL(IJ) > Zero)) then
        COMPARE = -OneHalf
      else if ((DV(IJ) >= OSTHR+dlt) .and. (DV(IJ) > Zero)) then
        COMPARE = -2.5_wp
      end if

      if (abs(COMPARE) >= TOLERANCE) then
        I_PRINT_HEADER = I_PRINT_HEADER+1
        if (I_PRINT_HEADER == 1) then
          write(u6,*)
          write(u6,*) ' Problematic transitions have been found'
          write(u6,*)
          write(u6,39) 'From','To','Difference (%)','Osc. st. (len.)','Osc. st. (vel.)'
          write(u6,40)
        end if
        if (COMPARE >= Zero) then
          write(u6,38) I,J,COMPARE*100.0_wp,DL(IJ),DV(IJ)
        else if (COMPARE >= -Two) then
          write(u6,36) I,J,DL(IJ),'below threshold'
        else
          write(u6,37) I,J,'below threshold',DV(IJ)
        end if
      end if

    end do ! L_
  end do ! K_
  if (I_PRINT_HEADER == 0) then
    write(u6,*)
    write(u6,*) 'No problematic oscillator strengths above the tolerance ',TOLERANCE,' have been found'
    write(u6,*)
  else
    write(u6,40)
    write(u6,*)
    write(u6,*) 'Number of problematic transitions = ',I_PRINT_HEADER
    write(u6,*)
  end if
  call CollapseOutput(0,'Length and velocity gauge comparison (spin-free states):')
  write(u6,*)
end if

! Free the memory

call mma_deallocate(DV)
call mma_deallocate(DL)

! CALCULATION OF THE QUADRUPOLE TRANSITION STRENGTHS

SECORD = 0

! We will first allocate a matrix for the total of the second order wave vector

call mma_allocate(TOT2K,NSTATE,NSTATE,Label='TOT2K')
TOT2K(:,:) = Zero

! Magnetic-Dipole - Magnetic-Dipole transitions

! Magnetic-Dipole

! DEBUG
!IPRDX_TEMP = IPRDX
!IPRDY_TEMP = IPRDY
!IPRDZ_TEMP = IPRDZ
! DEBUG END
IPRDX = 0
IPRDY = 0
IPRDZ = 0

IFANYD = 0
do IPROP=1,NPROP
  if (PNAME(IPROP) == 'ANGMOM') then
    IFANYD = 1
    if (ICOMP(IPROP) == 1) IPRDX = IPROP
    if (ICOMP(IPROP) == 2) IPRDY = IPROP
    if (ICOMP(IPROP) == 3) IPRDZ = IPROP
  end if
end do

if (IFANYD /= 0) then

  ! Only print the part calculated

  if (QIALL) then
    write(u6,*)
    call CollapseOutput(1,'Magnetic-Dipole - Magnetic-Dipole transition strengths (spin-free states):')
    write(u6,'(3X,A)') '--------------------------------------------------------------------------'
    if (OSTHR2 > Zero) then
      write(u6,30) 'for osc. strength at least',OSTHR2
      write(u6,*)
    end if
    write(u6,31) 'From','To','Osc. strength'
    write(u6,35)
  end if

  do K_=1,IEND
    I = IndexE(K_)
    do L_=JSTART,NSTATE
      J = IndexE(L_)
      EDIFF = ENERGY(J)-ENERGY(I)
      if (EDIFF > Zero) then

        DX2 = Zero
        DY2 = Zero
        DZ2 = Zero

        if (IPRDX > 0) DX2 = PROP(J,I,IPRDX)**2
        if (IPRDY > 0) DY2 = PROP(J,I,IPRDY)**2
        if (IPRDZ > 0) DZ2 = PROP(J,I,IPRDZ)**2

        F = (DX2+DY2+DZ2)*EDIFF*ONEOVER6C2
        ! Add it to the total
        TOT2K(J,I) = TOT2K(J,I)+F
        if (abs(F) >= OSTHR2) then
          !write(u6,*) ' value at distance'
          if (QIALL) write(u6,33) I,J,F
        end if

        ! Debug to move along z. Change DX and DY (1 Aangstrom)

        !do I=1,9
        !  AA = ZVAL(I)/Angstrom*EDIFF!/(Two*c_in_au)
        !  ! z-direction
        !  DX2 = (PROP(J,I,IPRDX)-AA*PROP(J,I,IPRDY_TEMP))**2
        !  DY2 = (PROP(J,I,IPRDY)+AA*PROP(J,I,IPRDX_TEMP))**2
        !  F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
        !  ! x-direction
        !  DZ2 = (PROP(J,I,IPRDZ)+AA*PROP(J,I,IPRDY_TEMP))**2
        !  DY2 = (PROP(J,I,IPRDY)-AA*PROP(J,I,IPRDZ_TEMP))**2
        !  F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
        !  ! y-direction
        !  DX2 = (PROP(J,I,IPRDX)+AA*PROP(J,I,IPRDZ_TEMP))**2
        !  DZ2 = (PROP(J,I,IPRDZ)-AA*PROP(J,I,IPRDX_TEMP))**2
        !  F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
        !  if (ABS(F) >= OSTHR2) then
        !    write(u6,*) ' moved value ',ZVAL(I)
        !    write(u6,'(5X,2I5,5X,G16.8)') I,J,F
        !  end if
        !end do
      end if
    end do
  end do
  if (QIALL) then
    write(u6,35)

    call CollapseOutput(0,'Magnetic-Dipole - Magnetic-Dipole transition strengths (spin-free states):')
  end if
  ! Magnetic-dipole - Magnetic-dipole calculated
  SECORD(1) = 1
end if

! Electric-Quadrupole Electric-Quadrupole transitions

IPRDXX = 0
IPRDXY = 0
IPRDXZ = 0
IPRDYY = 0
IPRDYZ = 0
IPRDZZ = 0

IFANYD = 0
do IPROP=1,NPROP
  if (PNAME(IPROP) == 'MLTPL  2') then
    IFANYD = 1
    if (ICOMP(IPROP) == 1) IPRDXX = IPROP
    if (ICOMP(IPROP) == 2) IPRDXY = IPROP
    if (ICOMP(IPROP) == 3) IPRDXZ = IPROP
    if (ICOMP(IPROP) == 4) IPRDYY = IPROP
    if (ICOMP(IPROP) == 5) IPRDYZ = IPROP
    if (ICOMP(IPROP) == 6) IPRDZZ = IPROP
  end if
end do

if (IFANYD /= 0) then
  if (QIALL) then
    write(u6,*)
    call CollapseOutput(1,'Quadrupole transition strengths (spin-free states):')
    write(u6,'(3X,A)') '---------------------------------------------------'
    if (OSTHR2 > Zero) then
      write(u6,30) 'for osc. strength at least',OSTHR2
      write(u6,*)
    end if
    write(u6,31) 'From','To','Osc. strength'
    write(u6,35)
  end if

  do K_=1,IEND
    I = IndexE(K_)
    do L_=JSTART,NSTATE
      J = IndexE(L_)
      EDIFF = ENERGY(J)-ENERGY(I)
      if (EDIFF > Zero) then

        EDIFF3 = EDIFF**3

        DXX = Zero
        DYY = Zero
        DZZ = Zero
        DXY = Zero
        DXZ = Zero
        DYZ = Zero
        if (IPRDXX > 0) DXX = PROP(J,I,IPRDXX)
        if (IPRDYY > 0) DYY = PROP(J,I,IPRDYY)
        if (IPRDZZ > 0) DZZ = PROP(J,I,IPRDZZ)
        if (IPRDXY > 0) DXY = PROP(J,I,IPRDXY)
        if (IPRDXZ > 0) DXZ = PROP(J,I,IPRDXZ)
        if (IPRDYZ > 0) DYZ = PROP(J,I,IPRDYZ)

        DXX2 = DXX**2
        DYY2 = DYY**2
        DZZ2 = DZZ**2
        FXX = ONEOVER30C*EDIFF3*(DXX2)
        FYY = ONEOVER30C*EDIFF3*(DYY2)
        FZZ = ONEOVER30C*EDIFF3*(DZZ2)

        DXY2 = DXY**2
        DXZ2 = DXZ**2
        DYZ2 = DYZ**2
        FXY = ONEOVER10C*EDIFF3*(DXY2)
        FXZ = ONEOVER10C*EDIFF3*(DXZ2)
        FYZ = ONEOVER10C*EDIFF3*(DYZ2)

        DXXDYY = DXX*DYY
        DXXDZZ = DXX*DZZ
        DYYDZZ = DYY*DZZ
        FXXFYY = -ONEOVER30C*EDIFF3*(DXXDYY)
        FXXFZZ = -ONEOVER30C*EDIFF3*(DXXDZZ)
        FYYFZZ = -ONEOVER30C*EDIFF3*(DYYDZZ)

        F = FXX+FXY+FXZ+FYY+FYZ+FZZ+FXXFYY+FXXFZZ+FYYFZZ
        ! Add it to the total
        TOT2K(J,I) = TOT2K(J,I)+F

        if (abs(F) >= OSTHR2) then
          if (QIALL) write(u6,33) I,J,F
        end if

        ! Debug to move along z. Change DZZ, DXZ and DYZ

        !do I=1,9
        !  DZZ2 = (PROP(J,I,IPRDZZ)+ZVAL(I)/Angstrom*Two*PROP(J,I,IPRDZ_TEMP))**2
        !  DXZ2 = (PROP(J,I,IPRDXZ)+ZVAL(I)/Angstrom*PROP(J,I,IPRDX_TEMP))**2
        !  DYZ2 = (PROP(J,I,IPRDYZ)+ZVAL(I)/Angstrom*PROP(J,I,IPRDY_TEMP))**2
        !  FZZ = ONEOVER30C*EDIFF3*(DZZ2)
        !  FXZ = ONEOVER10C*EDIFF3*(DXZ2)
        !  FYZ = ONEOVER10C*EDIFF3*(DYZ2)
        !  DXXDZZ = PROP(J,I,IPRDXX)*(PROP(J,I,IPRDZZ)+ZVAL(I)/Angstrom*Two*PROP(J,I,IPRDZ_TEMP))
        !  DYYDZZ = PROP(J,I,IPRDYY)*(PROP(J,I,IPRDZZ)+ZVAL(I)/Angstrom*Two*PROP(J,I,IPRDZ_TEMP))
        !  FXXFZZ = -ONEOVER30C*EDIFF3*(DXXDZZ)
        !  FYYFZZ = -ONEOVER30C*EDIFF3*(DYYDZZ)
        !  F = FXX+FXY+FXZ+FYY+FYZ+FZZ+FXXFYY+FXXFZZ+FYYFZZ
        !  if (ABS(F) >= OSTHR2) then
        !    write(u6,*) ' moved value ',ZVAL(I)
        !    write(u6,'(5X,2I5,5X,G16.8)') I,J,F
        !  end if
        !end do
      end if
    end do
  end do
  if (QIALL) then
    write(u6,35)

    call CollapseOutput(0,'Quadrupole transition strengths (spin-free states):')
  end if
  SECORD(2) = 1
end if

!Electric-Dipole Electric-Octupole transitions

! Octupole
! This is a real symmetric rank 3 tensor so only 10 and not 27 is needed
! The order which comes in

! DEBUG
!IPRDXX_TEMP = IPRDXX
!IPRDXY_TEMP = IPRDXY
!IPRDXZ_TEMP = IPRDXZ
!IPRDYY_TEMP = IPRDYY
!IPRDYZ_TEMP = IPRDYZ
!IPRDZZ_TEMP = IPRDZZ
! DEBUG END
IPRDXXX = 0
IPRDXXY = 0
IPRDXXZ = 0

!IPRDXYX = 0
!IPRDXYY = 0 ! YYX These are the same due to symmetry
!IPRDXYZ = 0 ! Not present

!IPRDXZX = 0
!IPRDXZY = 0
!IPRDXZZ = 0 ! ZZX

!IPRDYXX = 0
!IPRDYXY = 0
!IPRDYXZ = 0

IPRDYYX = 0 ! Taking the XYY order
IPRDYYY = 0 !
IPRDYYZ = 0 !

!IPRDYZX = 0
!IPRDYZY = 0
!IPRDYZZ = 0 ! ZZY

!IPRDZXX = 0
!IPRDZXY = 0
!IPRDZXZ = 0

!IPRDZYX = 0
!IPRDZYY = 0
!IPRDZYZ = 0

IPRDZZX = 0 ! Taking order from XZZ
IPRDZZY = 0 ! Taking order from YZZ
IPRDZZZ = 0
! Dipole
IPRDX = 0
IPRDY = 0
IPRDZ = 0

IFANYD = 0
do IPROP=1,NPROP
  if (PTYPE(IPROP)(5:8) /= 'SING') cycle
  if (PNAME(IPROP) == 'MLTPL  1') then
    if (ICOMP(IPROP) == 1) IPRDX = IPROP
    if (ICOMP(IPROP) == 2) IPRDY = IPROP
    if (ICOMP(IPROP) == 3) IPRDZ = IPROP
  else if (PNAME(IPROP) == 'MLTPL  3') then
    IFANYD = 1
    if (ICOMP(IPROP) == 1) IPRDXXX = IPROP
    if (ICOMP(IPROP) == 2) IPRDXXY = IPROP
    if (ICOMP(IPROP) == 3) IPRDXXZ = IPROP
    if (ICOMP(IPROP) == 4) IPRDYYX = IPROP ! Changed from XYY
    !if (ICOMP(IPROP) == 5) IPRDXYZ = IPROP
    if (ICOMP(IPROP) == 6) IPRDZZX = IPROP ! Changed from XZZ
    if (ICOMP(IPROP) == 7) IPRDYYY = IPROP
    if (ICOMP(IPROP) == 8) IPRDYYZ = IPROP
    if (ICOMP(IPROP) == 9) IPRDZZY = IPROP ! Changed from YZZ
    if (ICOMP(IPROP) == 10) IPRDZZZ = IPROP
  end if
end do
! Sanity check. Only check that dipole are there
! since it will give problems the other way when
! only calculating dipole transitions
if (((IPRDXXX > 0) .or. (IPRDYYX > 0) .or. (IPRDZZX > 0)) .and. (IPRDX <= 0)) then
  write(u6,*) ' Remember to include both Dipole and Octupole'
  call ABEND()
end if
if (((IPRDXXY > 0) .or. (IPRDYYY > 0) .or. (IPRDZZY > 0)) .and. (IPRDY <= 0)) then
  write(u6,*) ' Remember to include both Dipole and Octupole'
  call ABEND()
end if
if (((IPRDXXZ > 0) .or. (IPRDYYZ > 0) .or. (IPRDZZZ > 0)) .and. (IPRDZ <= 0)) then
  write(u6,*) ' Remember to include both Dipole and Octupole'
  call ABEND()
end if

if (IFANYD /= 0) then
  if (QIALL) then
    write(u6,*)
    call CollapseOutput(1,'Electric-Dipole - Electric-Octupole transition strengths (spin-free states):')
    write(u6,'(3X,A)') '----------------------------------------------------------------------------'
    if (OSTHR2 > Zero) then
      write(u6,30) 'for osc. strength at least',OSTHR2
      write(u6,*)
    end if
    write(u6,31) 'From','To','Osc. strength'
    write(u6,35)
  end if

  do K_=1,IEND
    I = IndexE(K_)
    do L_=JSTART,NSTATE
      J = IndexE(L_)
      EDIFF = ENERGY(J)-ENERGY(I)
      if (EDIFF > Zero) then

        EDIFF3 = EDIFF**3

        DXXXDX = Zero
        DYYXDX = Zero
        DZZXDX = Zero
        if (IPRDXXX > 0) DXXXDX = PROP(J,I,IPRDXXX)*PROP(J,I,IPRDX)
        if (IPRDYYX > 0) DYYXDX = PROP(J,I,IPRDYYX)*PROP(J,I,IPRDX)
        if (IPRDZZX > 0) DZZXDX = PROP(J,I,IPRDZZX)*PROP(J,I,IPRDX)
        FXXX = TWOOVERM45C*EDIFF3*(DXXXDX)
        FYYX = TWOOVERM45C*EDIFF3*(DYYXDX)
        FZZX = TWOOVERM45C*EDIFF3*(DZZXDX)

        DXXYDY = Zero
        DYYYDY = Zero
        DZZYDY = Zero
        if (IPRDXXY > 0) DXXYDY = PROP(J,I,IPRDXXY)*PROP(J,I,IPRDY)
        if (IPRDYYY > 0) DYYYDY = PROP(J,I,IPRDYYY)*PROP(J,I,IPRDY)
        if (IPRDZZY > 0) DZZYDY = PROP(J,I,IPRDZZY)*PROP(J,I,IPRDY)
        FXXY = TWOOVERM45C*EDIFF3*(DXXYDY)
        FYYY = TWOOVERM45C*EDIFF3*(DYYYDY)
        FZZY = TWOOVERM45C*EDIFF3*(DZZYDY)

        DXXZDZ = Zero
        DYYZDZ = Zero
        DZZZDZ = Zero
        if (IPRDXXZ > 0) DXXZDZ = PROP(J,I,IPRDXXZ)*PROP(J,I,IPRDZ)
        if (IPRDYYZ > 0) DYYZDZ = PROP(J,I,IPRDYYZ)*PROP(J,I,IPRDZ)
        if (IPRDZZZ > 0) DZZZDZ = PROP(J,I,IPRDZZZ)*PROP(J,I,IPRDZ)
        FXXZ = TWOOVERM45C*EDIFF3*(DXXZDZ)
        FYYZ = TWOOVERM45C*EDIFF3*(DYYZDZ)
        FZZZ = TWOOVERM45C*EDIFF3*(DZZZDZ)

        F = FXXX+FYYX+FZZX+FXXY+FYYY+FZZY+FXXZ+FYYZ+FZZZ
        ! Add it to the total
        TOT2K(J,I) = TOT2K(J,I)+F

        if (abs(F) >= OSTHR2) then
          if (QIALL) write(u6,33) I,J,F
        end if

        ! Debug to move along z. Change DZZX,DZZY,DXXZ,DYYZ and DZZZ

        !do I=1,9
        !  DZZXDX = (PROP(J,I,IPRDZZX)+ZVAL(I)/Angstrom*Two*PROP(J,I,IPRDXZ_TEMP)+(ZVAL(I)/Angstrom)**2*PROP(J,I,IPRDX_TEMP))* &
        !           PROP(J,I,IPRDX)
        !  FZZX = TWOOVERM45C*EDIFF3*(DZZXDX)
        !  DZZYDY = (PROP(J,I,IPRDZZY)+ZVAL(I)/Angstrom*Two*PROP(J,I,IPRDYZ_TEMP)+(ZVAL(I)/Angstrom)**2*PROP(J,I,IPRDY_TEMP))* &
        !           PROP(J,I,IPRDY)
        !  FZZY = TWOOVERM45C*EDIFF3*(DZZYDY)
        !  DXXZDZ = (PROP(J,I,IPRDXXZ)+ZVAL(I)/Angstrom*PROP(J,I,IPRDXX_TEMP))*PROP(J,I,IPRDZ)
        !  DYYZDZ = (PROP(J,I,IPRDYYZ)+ZVAL(I)/Angstrom*PROP(J,I,IPRDYY_TEMP))*PROP(J,I,IPRDZ)
        !  DZZZDZ = (PROP(J,I,IPRDZZZ)+ZVAL(I)/Angstrom*Three*PROP(J,I,IPRDZZ_TEMP)+ &
        !           (ZVAL(I)/Angstrom)**2*Three*PROP(J,I,IPRDZ_TEMP))*PROP(J,I,IPRDZ)
        !  FXXZ = TWOOVERM45C*EDIFF3*(DXXZDZ)
        !  FYYZ = TWOOVERM45C*EDIFF3*(DYYZDZ)
        !  FZZZ = TWOOVERM45C*EDIFF3*(DZZZDZ)
        !  F = FXXX+FYYX+FZZX+FXXY+FYYY+FZZY+FXXZ+FYYZ+FZZZ
        !  if (ABS(F) >= OSTHR2) then
        !    write(u6,*) ' moved value ',ZVAL(I)
        !    write(u6,'(5X,2I5,5X,G16.8)') I,J,F
        !  end if
        !end do
      end if
    end do
  end do
  if (QIALL) then
    write(u6,35)

    call CollapseOutput(0,'Electric-Dipole - Electric-Octupole transition strengths (spin-free states):')
  end if
  SECORD(3) = 1
end if

! Electric-Dipole - Magnetic-Quadrupole transitions

! Magnetic-Quadrupole
IPRDXX = 0
IPRDXY = 0
IPRDXZ = 0

IPRDYX = 0
IPRDYY = 0
IPRDYZ = 0

IPRDZX = 0
IPRDZY = 0
IPRDZZ = 0
! Electric-Dipole
IPRDX = 0
IPRDY = 0
IPRDZ = 0

IFANYD = 0
do IPROP=1,NPROP
  if (PTYPE(IPROP)(5:8) /= 'SING') cycle
  if (PNAME(IPROP) == 'MLTPL  1') then
    if (ICOMP(IPROP) == 1) IPRDX = IPROP
    if (ICOMP(IPROP) == 2) IPRDY = IPROP
    if (ICOMP(IPROP) == 3) IPRDZ = IPROP
  else if (PNAME(IPROP) == 'OMQ') then
    IFANYD = 1
    if (ICOMP(IPROP) == 1) IPRDXX = IPROP
    if (ICOMP(IPROP) == 2) IPRDXY = IPROP
    if (ICOMP(IPROP) == 3) IPRDXZ = IPROP

    if (ICOMP(IPROP) == 4) IPRDYX = IPROP
    if (ICOMP(IPROP) == 5) IPRDYY = IPROP
    if (ICOMP(IPROP) == 6) IPRDYZ = IPROP

    if (ICOMP(IPROP) == 7) IPRDZX = IPROP
    if (ICOMP(IPROP) == 8) IPRDZY = IPROP
    if (ICOMP(IPROP) == 9) IPRDZZ = IPROP
  end if
end do
! Sanity check. Only check that dipole are there
! since it will give problems the other way when
! only calculating dipole transitions
if (((IPRDYZ > 0) .or. (IPRDZY > 0)) .and. (IPRDX <= 0)) then
  write(u6,*) ' Remember to include both Dipole and Quadrupole'
  call ABEND()
end if
if (((IPRDZX > 0) .or. (IPRDXZ > 0)) .and. (IPRDY <= 0)) then
  write(u6,*) ' Remember to include both Dipole and Quadrupole'
  call ABEND()
end if
if (((IPRDXY > 0) .or. (IPRDYX > 0)) .and. (IPRDZ <= 0)) then
  write(u6,*) ' Remember to include both Dipole and Quadrupole'
  call ABEND()
end if

if (IFANYD /= 0) then
  if (QIALL) then
    write(u6,*)
    call CollapseOutput(1,'Electric-Dipole - Magnetic-Quadrupole transition strengths (spin-free states):')
    write(u6,'(3X,A)') '-----------------------------------------------------------------------------'

    if (OSTHR2 > Zero) then
      write(u6,30) 'for osc. strength at least',OSTHR2
      write(u6,*)
    end if
    write(u6,31) 'From','To','Osc. strength'
    write(u6,35)
  end if

  do K_=1,IEND
    I = IndexE(K_)
    do L_=JSTART,NSTATE
      J = IndexE(L_)
      EDIFF = ENERGY(J)-ENERGY(I)
      if (EDIFF > Zero) then

        EDIFF2 = EDIFF**2

        DXYDZ = Zero
        DYXDZ = Zero
        if (IPRDXY > 0) DXYDZ = PROP(J,I,IPRDXY)*PROP(J,I,IPRDZ)
        if (IPRDXY > 0) DYXDZ = PROP(J,I,IPRDYX)*PROP(J,I,IPRDZ)
        FXY = ONEOVER9C2*EDIFF2*(DXYDZ)
        FYX = -ONEOVER9C2*EDIFF2*(DYXDZ)

        DZXDY = Zero
        DXZDY = Zero
        if (IPRDZX > 0) DZXDY = PROP(J,I,IPRDZX)*PROP(J,I,IPRDY)
        if (IPRDXZ > 0) DXZDY = PROP(J,I,IPRDXZ)*PROP(J,I,IPRDY)
        FZX = ONEOVER9C2*EDIFF2*(DZXDY)
        FXZ = -ONEOVER9C2*EDIFF2*(DXZDY)

        DYZDX = Zero
        DZYDX = Zero
        if (IPRDYZ > 0) DYZDX = PROP(J,I,IPRDYZ)*PROP(J,I,IPRDX)
        if (IPRDZY > 0) DZYDX = PROP(J,I,IPRDZY)*PROP(J,I,IPRDX)
        FYZ = ONEOVER9C2*EDIFF2*(DYZDX)
        FZY = -ONEOVER9C2*EDIFF2*(DZYDX)

        F = FYX+FXY+FZX+FXZ+FYZ+FZY
        ! Add it to the total
        TOT2K(J,I) = TOT2K(J,I)+F

        if (abs(F) >= OSTHR2) then
          if (QIALL) write(u6,33) I,J,F
        end if

        ! Debug to move along z.

        !ONEOVER3C = One/(Three*c_in_au)
        !DYXDZ = (PROP(J,I,IPRDYX)*ONEOVER3C+PROP(J,I,IPRDXX_TEMP)*ONEOVER3C*EDIFF/Angstrom)*PROP(J,I,IPRDZ)
        !DXYDZ = (PROP(J,I,IPRDXY)*ONEOVER3C-PROP(J,I,IPRDYY_TEMP)*ONEOVER3C*EDIFF/Angstrom)*PROP(J,I,IPRDZ)
        !write(u6,*) 'YX,XY moved',PROP(J,I,IPRDYX)+PROP(J,I,IPRDXX_TEMP)*EDIFF/Angstrom, &
        !            PROP(J,I,IPRDXY)-PROP(J,I,IPRDYY_TEMP)*EDIFF/Angstrom
        !FXY = ONEOVER3C*EDIFF2*(DXYDZ)
        !FYX = -ONEOVER3C*EDIFF2*(DYXDZ)
        !DZXDY = PROP(J,I,IPRDZX)*ONEOVER3C*PROP(J,I,IPRDY) !independent
        !DXZDY = (PROP(J,I,IPRDXZ)*ONEOVER3C-PROP(J,I,IPRDYZ_TEMP)*ONEOVER3C*EDIFF/Angstrom+ &
        !         PROP(J,I,IPRDY_TEMP)*Two*ONEOVER3C*EDIFF/Angstrom**2)*PROP(J,I,IPRDY) ! skipped magnetic dipole
        !write(u6,*) 'ZX,XZ moved',PROP(J,I,IPRDZX), &
        !            PROP(J,I,IPRDXZ)-PROP(J,I,IPRDYZ_TEMP)*EDIFF/Angstrom+PROP(J,I,IPRDY_TEMP)*Two*EDIFF/Angstrom**2
        !FZX = ONEOVER3C*EDIFF2*(DZXDY)
        !FXZ = -ONEOVER3C*EDIFF2*(DXZDY)
        !DYZDX = (PROP(J,I,IPRDYZ)*ONEOVER3C+PROP(J,I,IPRDXZ_TEMP)*ONEOVER3C*EDIFF/Angstrom- &
        !         PROP(J,I,IPRDX_TEMP)*Two*ONEOVER3C*EDIFF/Angstrom**2)*PROP(J,I,IPRDX)
        !DZYDX = PROP(J,I,IPRDZY)*ONEOVER3C*PROP(J,I,IPRDX)
        !write(u6,*) 'YZ,ZY moved',PROP(J,I,IPRDYZ)+ &
        !             PROP(J,I,IPRDXZ_TEMP)*EDIFF/Angstrom-PROP(J,I,IPRDX_TEMP)*Two*EDIFF/Angstrom**2,PROP(J,I,IPRDZY)
        !FYZ = ONEOVER3C*EDIFF2*(DYZDX)
        !FZY = -ONEOVER3C*EDIFF2*(DZYDX)
        ! The new diagonal ones?
        !F = FYX+FXY+FZX+FXZ+FYZ+FZY
        !if (ABS(F) >= OSTHR2) then
        !  write(u6,*) ' The moved value'
        !  write(u6,'(5X,2I5,5X,G16.8)') I,J,F
        !end if
        ! End debug
      end if
    end do
  end do
  if (QIALL) then
    write(u6,35)

    call CollapseOutput(0,'Electric-Dipole - Magnetic-Quadrupole transition strengths (spin-free states):')
  end if
  SECORD(4) = 1
end if

! Now write out the total

! Add it to the total

I2TOT = count(SECORD(:) == 1)
if (I2TOT >= 1) then
  if (SECORD(1) == 0) write(u6,*) 'Magnetic-dipole - magnetic-dipole not included'
  if (SECORD(2) == 0) write(u6,*) 'Electric-quadrupole - electric-quadrupole not included'
  if (SECORD(3) == 0) write(u6,*) 'Electric-dipole - electric-octupole not included'
  if (SECORD(4) == 0) write(u6,*) 'Electric-dipole - magnetic-quadrupole not included'
  iPrint = 0
  do K_=1,IEND
    I = IndexE(K_)
    do L_=JSTART,NSTATE
      J = IndexE(L_)
      EDIFF = ENERGY(J)-ENERGY(I)
      if (EDIFF > Zero) then

        F = TOT2K(J,I)
        if (abs(F) >= OSTHR2) then
          if (iPrint == 0) then
            write(u6,*)
            call CollapseOutput(1,'Second-order contribution to the transition strengths (spin-free states):')
            write(u6,'(3X,A)') '-------------------------------------------------------------------------'

            if (OSTHR2 > Zero) then
              write(u6,30) 'for osc. strength at least',OSTHR2
              write(u6,*)
            end if
            write(u6,31) 'From','To','Osc. strength'
            write(u6,35)
            iPrint = 1
          end if
          write(u6,33) I,J,F
          call Add_Info('TMS(SF,2nd)',[F],1,6)
        end if
      end if
    end do
  end do
  if (iPrint == 1) then
    write(u6,35)
    call CollapseOutput(0,'Second-order contribution to the transition strengths (spin-free states):')
  end if
end if
! release the memory again
call mma_deallocate(TOT2K)

if (DOCD) then
  ! Lasse 2019
  ! New CD here with electric dipole and magnetic-dipole - velocity gauge
  IPRDXD = 0
  IPRDYD = 0
  IPRDZD = 0
  IPRDXM = 0
  IPRDYM = 0
  IPRDZM = 0
  IPRQXX = 0
  IPRQXY = 0
  IPRQXZ = 0
  IPRQYY = 0
  IPRQYZ = 0
  IPRQZZ = 0

  IFANYD = 0
  IFANYM = 0
  IFANYQ = 0
  do IPROP=1,NPROP
    if (PNAME(IPROP) == 'VELOCITY') then
      IFANYD = 1
      if (ICOMP(IPROP) == 1) IPRDXD = IPROP
      if (ICOMP(IPROP) == 2) IPRDYD = IPROP
      if (ICOMP(IPROP) == 3) IPRDZD = IPROP
    else if (PNAME(IPROP) == 'ANGMOM') then
      IFANYM = 1
      if (ICOMP(IPROP) == 1) IPRDXM = IPROP
      if (ICOMP(IPROP) == 2) IPRDYM = IPROP
      if (ICOMP(IPROP) == 3) IPRDZM = IPROP
    else if (PNAME(IPROP) == 'MLTPV  2') then
      IFANYQ = 1
      if (ICOMP(IPROP) == 1) IPRQXX = IPROP
      if (ICOMP(IPROP) == 2) IPRQXY = IPROP
      if (ICOMP(IPROP) == 3) IPRQXZ = IPROP
      if (ICOMP(IPROP) == 4) IPRQYY = IPROP
      if (ICOMP(IPROP) == 5) IPRQYZ = IPROP
      if (ICOMP(IPROP) == 6) IPRQZZ = IPROP
    end if
  end do

  if ((IFANYD /= 0) .and. (IFANYM /= 0)) then

    ! Only print the part calculated

    write(u6,*)
    call CollapseOutput(1,'Circular Dichroism - velocity gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (spin-free '// &
                        'states):')
    write(u6,'(3X,A)') '----------------------------------------------------------------------------------------------------'// &
                       '--------'
    if (DO_SK) then
      write(u6,30) 'For red. rot. strength at least',RSTHR
    else
      write(u6,30) 'For isotropic red. rot. strength at least',RSTHR
    end if

    write(u6,*)

    if (Do_SK .and. (IFANYQ /= 0)) then
      nVec = nk_Vector
    else
      nVec = 1
    end if

    do iVec=1,nVec

      if (Do_SK .and. (IFANYQ /= 0)) then
        write(u6,*)
        write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3)
        write(u6,*)
        write(u6,31) 'From','To','Red. rot. str.'
      else
        write(u6,31) 'From','To','Red. rot. str.'
        if (IFANYQ /= 0) write(u6,44) 'Rxx','Rxy','Rxz','Ryy','Ryz','Rzz'
      end if
      write(u6,35)

      do K_=1,IEND
        I = IndexE(K_)
        do L_=JSTART,NSTATE
          J = IndexE(L_)
          EDIFF = ENERGY(J)-ENERGY(I)
          if (EDIFF > Zero) then

            ! R = e^2*hbar/(2*m^2*E) <J|p|I>.<I|l|J>
            !   = e^2*hbar/(2*m^2*E) -i*hbar*<J|nabla|I>.-i*hbar*<I|r x nabla|J>
            !   = e^2*hbar^3/(2*m^2*E) <J|nabla|I>.<J|r x nabla|I>

            RXX = Zero
            RYY = Zero
            RZZ = Zero
            if ((IPRDXD > 0) .and. (IPRDXM > 0)) RXX = PROP(J,I,IPRDXD)*PROP(J,I,IPRDXM)
            if ((IPRDYD > 0) .and. (IPRDYM > 0)) RYY = PROP(J,I,IPRDYD)*PROP(J,I,IPRDYM)
            if ((IPRDZD > 0) .and. (IPRDZM > 0)) RZZ = PROP(J,I,IPRDZD)*PROP(J,I,IPRDZM)
            R = RXX+RYY+RZZ
            R = R*Half/EDIFF*AU2REDR

            ! Compute full rotatory strength tensor
            ! (see Hansen and Bak, 10.1021/jp001899+)

            if (IFANYQ /= 0) then
              RXY = Zero
              RXZ = Zero
              RYX = Zero
              RYZ = Zero
              RZX = Zero
              RZY = Zero
              RXXY = Zero
              RXXZ = Zero
              RXYX = Zero
              RXYZ = Zero
              RXZX = Zero
              RXZY = Zero
              RXYY = Zero
              RYYX = Zero
              RYYZ = Zero
              RYZX = Zero
              RYZY = Zero
              RXZZ = Zero
              RYZZ = Zero
              RZZX = Zero
              RZZY = Zero
              if ((IPRDXD > 0) .and. (IPRDYM > 0)) RXY = PROP(J,I,IPRDXD)*PROP(J,I,IPRDYM)
              if ((IPRDXD > 0) .and. (IPRDZM > 0)) RXZ = PROP(J,I,IPRDXD)*PROP(J,I,IPRDZM)
              if ((IPRDYD > 0) .and. (IPRDXM > 0)) RYX = PROP(J,I,IPRDYD)*PROP(J,I,IPRDXM)
              if ((IPRDYD > 0) .and. (IPRDZM > 0)) RYZ = PROP(J,I,IPRDYD)*PROP(J,I,IPRDZM)
              if ((IPRDZD > 0) .and. (IPRDXM > 0)) RZX = PROP(J,I,IPRDZD)*PROP(J,I,IPRDXM)
              if ((IPRDZD > 0) .and. (IPRDYM > 0)) RZY = PROP(J,I,IPRDZD)*PROP(J,I,IPRDYM)
              if ((IPRQXX > 0) .and. (IPRDYD > 0)) RXXY = PROP(J,I,IPRQXX)*PROP(J,I,IPRDYD)
              if ((IPRQXX > 0) .and. (IPRDZD > 0)) RXXZ = PROP(J,I,IPRQXX)*PROP(J,I,IPRDZD)
              if ((IPRQXY > 0) .and. (IPRDXD > 0)) RXYX = PROP(J,I,IPRQXY)*PROP(J,I,IPRDXD)
              if ((IPRQXY > 0) .and. (IPRDZD > 0)) RXYZ = PROP(J,I,IPRQXY)*PROP(J,I,IPRDZD)
              if ((IPRQXZ > 0) .and. (IPRDXD > 0)) RXZX = PROP(J,I,IPRQXZ)*PROP(J,I,IPRDXD)
              if ((IPRQXZ > 0) .and. (IPRDYD > 0)) RXZY = PROP(J,I,IPRQXZ)*PROP(J,I,IPRDYD)
              if ((IPRQXY > 0) .and. (IPRDYD > 0)) RXYY = PROP(J,I,IPRQXY)*PROP(J,I,IPRDYD)
              if ((IPRQYY > 0) .and. (IPRDXD > 0)) RYYX = PROP(J,I,IPRQYY)*PROP(J,I,IPRDXD)
              if ((IPRQYY > 0) .and. (IPRDZD > 0)) RYYZ = PROP(J,I,IPRQYY)*PROP(J,I,IPRDZD)
              if ((IPRQYZ > 0) .and. (IPRDXD > 0)) RYZX = PROP(J,I,IPRQYZ)*PROP(J,I,IPRDXD)
              if ((IPRQYZ > 0) .and. (IPRDYD > 0)) RYZY = PROP(J,I,IPRQYZ)*PROP(J,I,IPRDYD)
              if ((IPRQXZ > 0) .and. (IPRDZD > 0)) RXZZ = PROP(J,I,IPRQXZ)*PROP(J,I,IPRDZD)
              if ((IPRQYZ > 0) .and. (IPRDZD > 0)) RYZZ = PROP(J,I,IPRQYZ)*PROP(J,I,IPRDZD)
              if ((IPRQZZ > 0) .and. (IPRDXD > 0)) RZZX = PROP(J,I,IPRQZZ)*PROP(J,I,IPRDXD)
              if ((IPRQZZ > 0) .and. (IPRDYD > 0)) RZZY = PROP(J,I,IPRQZZ)*PROP(J,I,IPRDYD)
              ! xx, xy, xz, yy, yz, zz
              Rtensor(1) = 0.75_wp*(RYY+RZZ+(RXYZ-RXZY))
              Rtensor(2) = -0.375_wp*(RXY+RYX+(RXXZ+RYZY-RXZX-RYYZ))
              Rtensor(3) = -0.375_wp*(RXZ+RZX+(RXYX+RZZY-RXXY-RYZZ))
              Rtensor(4) = 0.75_wp*(RXX+RZZ+(RYZX-RXYZ))
              Rtensor(5) = -0.375_wp*(RYZ+RZY+(RYYX+RXZZ-RXYY-RZZX))
              Rtensor(6) = 0.75_wp*(RXX+RYY+(RXZY-RYZX))
              call DSCAL_(6,AU2REDR/EDIFF,Rtensor,1)
              if (Do_SK) then
                ! k^T R k
                R = k_vector(1,iVec)**2*Rtensor(1)+k_vector(2,iVec)**2*Rtensor(4)+k_vector(3,iVec)**2*Rtensor(6)+ &
                    Two*k_vector(1,iVec)*k_vector(2,iVec)*Rtensor(2)+Two*k_vector(1,iVec)*k_vector(3,iVec)*Rtensor(3)+ &
                    Two*k_vector(2,iVec)*k_vector(3,iVec)*Rtensor(5)
              else
                write(u6,43) 'tensor: ',Rtensor(:)
              end if
            end if

            if (abs(R) > RSTHR) write(u6,33) I,J,R

            call Add_Info('CD_V(SF)',[R],1,6)
          end if
        end do
      end do
      write(u6,35)
    end do

    call CollapseOutput(0,'Circular Dichroism - velocity gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (spin-free '// &
                        'states):')
  end if

  ! Lasse 2019
  ! New CD here with electric dipole and magnetic-dipole - mixed gauge
  ! Usually refered to as the length gauge
  IPRDXD = 0
  IPRDYD = 0
  IPRDZD = 0
  IPRDXM = 0
  IPRDYM = 0
  IPRDZM = 0
  IPRQXX = 0
  IPRQXY = 0
  IPRQXZ = 0
  IPRQYY = 0
  IPRQYZ = 0
  IPRQZZ = 0

  IFANYD = 0
  IFANYM = 0
  IFANYQ = 0
  do IPROP=1,NPROP
    if (PTYPE(IPROP)(5:8) /= 'SING') cycle
    if (PNAME(IPROP) == 'MLTPL  1') then
      IFANYD = 1
      if (ICOMP(IPROP) == 1) IPRDXD = IPROP
      if (ICOMP(IPROP) == 2) IPRDYD = IPROP
      if (ICOMP(IPROP) == 3) IPRDZD = IPROP
    else if (PNAME(IPROP) == 'ANGMOM') then
      IFANYM = 1
      if (ICOMP(IPROP) == 1) IPRDXM = IPROP
      if (ICOMP(IPROP) == 2) IPRDYM = IPROP
      if (ICOMP(IPROP) == 3) IPRDZM = IPROP
    else if (PNAME(IPROP) == 'MLTPL  2') then
      IFANYQ = 1
      if (ICOMP(IPROP) == 1) IPRQXX = IPROP
      if (ICOMP(IPROP) == 2) IPRQXY = IPROP
      if (ICOMP(IPROP) == 3) IPRQXZ = IPROP
      if (ICOMP(IPROP) == 4) IPRQYY = IPROP
      if (ICOMP(IPROP) == 5) IPRQYZ = IPROP
      if (ICOMP(IPROP) == 6) IPRQZZ = IPROP
    end if
  end do

  if ((IFANYD /= 0) .and. (IFANYM /= 0)) then

    ! Only print the part calculated

    write(u6,*)
    call CollapseOutput(1,'Circular Dichroism - mixed gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (spin-free '// &
                        'states):')
    write(u6,'(3X,A)') '---------------------------------------------------------------------------------------------------------'
    write(u6,*)
    write(u6,*) ' WARNING WARNING WARNING !!!'
    write(u6,*)
    write(u6,*) ' Circular Dichroism in the mixed gauge'
    write(u6,*) ' is NOT origin independent - check your results'
    if (DO_SK) then
      write(u6,30) 'For red. rot. strength at least',RSTHR
    else
      write(u6,30) 'For isotropic red. rot. strength at least',RSTHR
    end if
    write(u6,*)

    if (Do_SK .and. (IFANYQ /= 0)) then
      nVec = nk_Vector
    else
      nVec = 1
    end if

    do iVec=1,nVec

      if (Do_SK .and. (IFANYQ /= 0)) then
        write(u6,*)
        write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3)
        write(u6,*)
        write(u6,31) 'From','To','Red. rot. str.'
      else
        write(u6,31) 'From','To','Red. rot. str.'
        if (IFANYQ /= 0) write(u6,44) 'Rxx','Rxy','Rxz','Ryy','Ryz','Rzz'
      end if
      write(u6,35)

      do K_=1,IEND
        I = IndexE(K_)
        do L_=JSTART,NSTATE
          J = IndexE(L_)
          EDIFF = ENERGY(J)-ENERGY(I)
          if (EDIFF > Zero) then

            ! R = -i*e^2/(2*m) <J|r|I>.<I|l|J>
            !   = -i*e^2/(2*m) <J|r|I>.-i*hbar*<I|r x nabla|J>
            !   = e^2*hbar/(2*m) <J|r|I>.<J|r x nabla|I>

            RXX = Zero
            RYY = Zero
            RZZ = Zero
            if ((IPRDXD > 0) .and. (IPRDXM > 0)) RXX = PROP(J,I,IPRDXD)*PROP(J,I,IPRDXM)
            if ((IPRDYD > 0) .and. (IPRDYM > 0)) RYY = PROP(J,I,IPRDYD)*PROP(J,I,IPRDYM)
            if ((IPRDZD > 0) .and. (IPRDZM > 0)) RZZ = PROP(J,I,IPRDZD)*PROP(J,I,IPRDZM)
            R = RXX+RYY+RZZ
            R = R*Half*AU2REDR

            ! Compute full rotatory strength tensor
            ! (see Hansen and Bak, 10.1021/jp001899+)

            if (IFANYQ /= 0) then
              RXY = Zero
              RXZ = Zero
              RYX = Zero
              RYZ = Zero
              RZX = Zero
              RZY = Zero
              RXXY = Zero
              RXXZ = Zero
              RXYX = Zero
              RXYZ = Zero
              RXZX = Zero
              RXZY = Zero
              RXYY = Zero
              RYYX = Zero
              RYYZ = Zero
              RYZX = Zero
              RYZY = Zero
              RXZZ = Zero
              RYZZ = Zero
              RZZX = Zero
              RZZY = Zero
              if ((IPRDXD > 0) .and. (IPRDYM > 0)) RXY = PROP(J,I,IPRDXD)*PROP(J,I,IPRDYM)
              if ((IPRDXD > 0) .and. (IPRDZM > 0)) RXZ = PROP(J,I,IPRDXD)*PROP(J,I,IPRDZM)
              if ((IPRDYD > 0) .and. (IPRDXM > 0)) RYX = PROP(J,I,IPRDYD)*PROP(J,I,IPRDXM)
              if ((IPRDYD > 0) .and. (IPRDZM > 0)) RYZ = PROP(J,I,IPRDYD)*PROP(J,I,IPRDZM)
              if ((IPRDZD > 0) .and. (IPRDXM > 0)) RZX = PROP(J,I,IPRDZD)*PROP(J,I,IPRDXM)
              if ((IPRDZD > 0) .and. (IPRDYM > 0)) RZY = PROP(J,I,IPRDZD)*PROP(J,I,IPRDYM)
              if ((IPRQXX > 0) .and. (IPRDYD > 0)) RXXY = PROP(J,I,IPRQXX)*PROP(J,I,IPRDYD)
              if ((IPRQXX > 0) .and. (IPRDZD > 0)) RXXZ = PROP(J,I,IPRQXX)*PROP(J,I,IPRDZD)
              if ((IPRQXY > 0) .and. (IPRDXD > 0)) RXYX = PROP(J,I,IPRQXY)*PROP(J,I,IPRDXD)
              if ((IPRQXY > 0) .and. (IPRDZD > 0)) RXYZ = PROP(J,I,IPRQXY)*PROP(J,I,IPRDZD)
              if ((IPRQXZ > 0) .and. (IPRDXD > 0)) RXZX = PROP(J,I,IPRQXZ)*PROP(J,I,IPRDXD)
              if ((IPRQXZ > 0) .and. (IPRDYD > 0)) RXZY = PROP(J,I,IPRQXZ)*PROP(J,I,IPRDYD)
              if ((IPRQXY > 0) .and. (IPRDYD > 0)) RXYY = PROP(J,I,IPRQXY)*PROP(J,I,IPRDYD)
              if ((IPRQYY > 0) .and. (IPRDXD > 0)) RYYX = PROP(J,I,IPRQYY)*PROP(J,I,IPRDXD)
              if ((IPRQYY > 0) .and. (IPRDZD > 0)) RYYZ = PROP(J,I,IPRQYY)*PROP(J,I,IPRDZD)
              if ((IPRQYZ > 0) .and. (IPRDXD > 0)) RYZX = PROP(J,I,IPRQYZ)*PROP(J,I,IPRDXD)
              if ((IPRQYZ > 0) .and. (IPRDYD > 0)) RYZY = PROP(J,I,IPRQYZ)*PROP(J,I,IPRDYD)
              if ((IPRQXZ > 0) .and. (IPRDZD > 0)) RXZZ = PROP(J,I,IPRQXZ)*PROP(J,I,IPRDZD)
              if ((IPRQYZ > 0) .and. (IPRDZD > 0)) RYZZ = PROP(J,I,IPRQYZ)*PROP(J,I,IPRDZD)
              if ((IPRQZZ > 0) .and. (IPRDXD > 0)) RZZX = PROP(J,I,IPRQZZ)*PROP(J,I,IPRDXD)
              if ((IPRQZZ > 0) .and. (IPRDYD > 0)) RZZY = PROP(J,I,IPRQZZ)*PROP(J,I,IPRDYD)
              ! xx, xy, xz, yy, yz, zz
              Rtensor(1) = 0.75_wp*(RYY+RZZ+EDIFF*(RXYZ-RXZY))
              Rtensor(2) = -0.375_wp*(RXY+RYX+EDIFF*(RXXZ+RYZY-RXZX-RYYZ))
              Rtensor(3) = -0.375_wp*(RXZ+RZX+EDIFF*(RXYX+RZZY-RXXY-RYZZ))
              Rtensor(4) = 0.75_wp*(RXX+RZZ+EDIFF*(RYZX-RXYZ))
              Rtensor(5) = -0.375_wp*(RYZ+RZY+EDIFF*(RYYX+RXZZ-RXYY-RZZX))
              Rtensor(6) = 0.75_wp*(RXX+RYY+EDIFF*(RXZY-RYZX))
              call DSCAL_(6,AU2REDR,Rtensor,1)
              if (Do_SK) then
                ! k^T R k
                R = k_vector(1,iVec)**2*Rtensor(1)+k_vector(2,iVec)**2*Rtensor(4)+k_vector(3,iVec)**2*Rtensor(6)+ &
                    Two*k_vector(1,iVec)*k_vector(2,iVec)*Rtensor(2)+Two*k_vector(1,iVec)*k_vector(3,iVec)*Rtensor(3)+ &
                    Two*k_vector(2,iVec)*k_vector(3,iVec)*Rtensor(5)
              else
                write(u6,43) 'tensor: ',Rtensor(:)
              end if
            end if

            if (abs(R) > RSTHR) write(u6,33) I,J,R

            call Add_Info('CD_M(SF)',[R],1,6)
          end if
        end do
      end do
      write(u6,35)
    end do

    call CollapseOutput(0,'Circular Dichroism - mixed gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (spin-free '// &
                        'states):')
  end if
end if
! CD end

! +++ J. Norell 12/7 - 2018
! Dyson amplitudes for (1-electron) ionization transitions
!+++ Bruno Tenorio, 2020. Added Corrected Dyson norms
! according to dysnorm subroutine.
if (DYSO) then
  DYSTHR = 1.0e-5_wp
  write(u6,*)
  call CollapseOutput(1,'Dyson amplitudes Biorth. corrected (spin-free states):')
  write(u6,'(3X,A)') '------------------------------------------------------'
  if (DYSTHR > Zero) then
    write(u6,30) 'for Dyson intensities at least',DYSTHR
    write(u6,30)
  end if
  write(u6,*) '       From      To        BE (eV)           Dyson intensity'
  write(u6,32)
  FMAX = Zero
  do I_=1,NSTATE
    I = IndexE(I_)
    do J_=1,NSTATE
      J = IndexE(J_)
      F = DYSAMPS2(I,J)*DYSAMPS2(I,J)
      EDIFF = auToeV*(ENERGY(J)-ENERGY(I))
      if (F > 1.0e-36_wp) then
        if (EDIFF > Zero) write(u6,'(A,I8,I8,F15.3,ES22.5)') '    ',I,J,EDIFF,F
      end if
    end do ! J
  end do ! I
  call CollapseOutput(0,'Dyson amplitudes Biorth. corrected (spin-free states):')
  write(u6,*)
end if
! +++ J. Norell

if (.not. Do_TMOM) then
  call FinishUp()
  return
end if

!***********************************************************************
!                                                                      *
!     Start of section for transition moments                          *
!                                                                      *
!***********************************************************************

! Find the section of transition moments in the property list.

! The operator is split in 4 different component, each with three
! elements corresponding to differentiation in the x, y, and z
! direction.

!***********************************************************************
!                                                                      *
!     Computation of the isotropic oscillator strength.                *
!                                                                      *
!***********************************************************************

#define _TIME_TMOM_
#ifdef _TIME_TMOM_
call CWTime(TCpu1,TWall1)
#endif
call mma_Allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
call mma_Allocate(TSDMZZ,nTDMZZ,Label='TSDMZZ')
call mma_Allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
nSCR = (NBST*(NBST+1))/2
call mma_allocate(SCR,nSCR,4,LABEL='SCR')

! Here we will use a Lebedev grid to integrate over all possible
! directions of the wave vector, k. The property integrals will be
! computed on the fly and traced with the density to generate the
! corresponding values in the PROP matrix.

! Find the slot on the one-electron file where we will store the
! on-the-fly generated property integrals.

IPRTMOM(:) = -1
do IPROP=1,NPROP
  if (PNAME(IPROP) == 'TMOM  RS') then
    if (IPRTMOM(0+ICOMP(IPROP)) == -1) IPRTMOM(0+ICOMP(IPROP)) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM  IS') then
    if (IPRTMOM(3+ICOMP(IPROP)) == -1) IPRTMOM(3+ICOMP(IPROP)) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM  RA') then
    if (IPRTMOM(6+ICOMP(IPROP)) == -1) IPRTMOM(6+ICOMP(IPROP)) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM  IA') then
    if (IPRTMOM(9+ICOMP(IPROP)) == -1) IPRTMOM(9+ICOMP(IPROP)) = IPROP
  end if
end do
if (any(IPRTMOM == -1)) return

! Initiate the Seward environment

nDiff = 0
call IniSew(.false.,nDiff)

! Generate the quadrature points.

if (Do_SK) then
  nQuad = 1
  call mma_Allocate(Rquad,4,nQuad,label='SK')
  nVec = nk_Vector
else
  call unitmat(Pax,3)
  ! In the spin-free case, oscillator and rotatory strengths for k and -k
  ! are equal, so we compute only half the quadrature points
  call Do_Lebedev(L_Eff,nQuad,Rquad,4)
  nVec = 1
end if
if (Do_Pol) call mma_allocate(pol_Vector,3,nVec*nQuad,Label='POL')

! Scratch for one-electron integrals

NIP = 4+(NBST*(NBST+1))/2
call mma_allocate(IP,NIP,Label='IP')
#ifdef _HDF5_

! Allocate vector to store all individual transition moments.
! We do this for
! all unique pairs I-J, I=/=J (NSTATE*(NSTATE-1)/2)
!     all k-vectors (nQuad or nVec)
!         we store:
!             the weight (1)
!             the k-vector (3)
!             the projected transition vector (real and imaginary parts) (2*3)

nIJ = nState*(nState-1)/2
ip_w = 1
ip_kvector = ip_w+1
ip_TMR = ip_kvector+3
ip_TMI = ip_TMR+3
nData = ip_TMI+3-1
call mma_allocate(Storage,nData,nQuad,nIJ,nVec,label='Storage')
call dCopy_(size(Storage),[Zero],0,Storage,1)
#endif
! MGD create the groups of indices
! Only with reduce loop to make things easier
TMOgroup = .false.
ngroup1 = IEND
ngroup2 = NSTATE-JSTART+1
nmax2 = 1
if (REDUCELOOP .and. (TMGr_thrs >= Zero)) then
  TMOgroup = .true.
  THRS = TMGr_thrs
  i = IndexE(IEND)
  RefEne = 0
  TAU = -1
  ngroup2 = 1
  do j_=JSTART,NSTATE
    j = IndexE(j_)
    if (ENERGY(J)-Refene > TAU) then
      NGROUP2 = NGROUP2+1
      Refene = ENERGY(J)
      ediff = Refene-ENERGY(I)
      TAU = ediff*THRS
    end if
  end do
  call mma_Allocate(TMOgrp2,NGROUP2,Label='TMOgrp2')
  ngroup2 = 0
  TAU = -1
  RefEne = 0
  do j_=JSTART,NSTATE
    j = IndexE(j_)
    if (ENERGY(J)-Refene > TAU) then
      NGROUP2 = NGROUP2+1
      TMOgrp2(NGROUP2) = J_
      Refene = ENERGY(J)
      ediff = Refene-ENERGY(I)
      TAU = ediff*THRS
    end if
  end do
  TMOgrp2(ngroup2+1) = NSTATE+1

  j = IndexE(JSTART)
  Refene = ENERGY(j)
  TAU = -1
  ngroup1 = 1
  do i_=IEND,1,-1
    i = IndexE(i_)
    if (Refene-ENERGY(i) > TAU) then
      ngroup1 = ngroup1+1
      Refene = energy(i)
      ediff = energy(j)-Refene
      Tau = ediff*THRS
    end if
  end do
  call mma_Allocate(TMOgrp1,NGROUP1,Label='TMOgrp1')
  Ntmp = Ngroup1
  Ngroup1 = Ngroup1-1
  Refene = ENERGY(j)
  TAU = -1
  do i_=IEND,1,-1
    i = IndexE(i_)
    if (Refene-ENERGY(i) > TAU) then
      TMOgrp1(ntmp) = i_+1
      ntmp = ntmp-1
      Refene = energy(i)
      ediff = energy(j)-Refene
      Tau = ediff*THRS
    end if
  end do
  TMOgrp1(1) = 1
  !write(u6,*) (TMOgrp1(i),i=1,ngroup1+1)
  !write(u6,*) (TMOgrp2(i),i=1,ngroup2+1)
  maxgrp1 = 0
  do i=1,ngroup1
    maxgrp1 = max(maxgrp1,TMOgrp1(i+1)-TMOgrp1(i))
  end do
  maxgrp2 = 0
  do i=1,ngroup2
    maxgrp2 = max(maxgrp2,TMOgrp2(i+1)-TMOgrp2(i))
  end do
  nmax2 = maxgrp1*maxgrp2
end if

! Array for printing contributions from different directions

call mma_allocate(RAW,NQUAD,6,nmax2,Label='RAW')
call mma_allocate(OSCSTR,2,nmax2,Label='OscStr')
call mma_allocate(Aux,8,nmax2,Label='Aux')

do iVec=1,nVec
  if (Do_SK) then
    Rquad(1:3,1) = k_Vector(:,iVec)
    Rquad(4,1) = One   ! Dummy weight
  end if

  iPrint = 0
  IJSO = 0
  do igrp=1,ngroup1
    do jgrp=1,ngroup2

      if (TMOgroup) then
        istart_ = TMOgrp1(igrp)
        iend_ = TMOgrp1(igrp+1)-1
        jstart_ = TMOgrp2(jgrp)
        jend_ = TMOgrp2(jgrp+1)-1
        EDIFF_ = (ENERGY(IndexE(jstart_))+ENERGY(IndexE(jend_))-ENERGY(IndexE(istart_))-ENERGY(IndexE(iend_)))*Half
      else
        istart_ = igrp
        iend_ = igrp
        jstart_ = jgrp+jstart-1
        jend_ = jgrp+jstart-1
        EDIFF_ = ENERGY(IndexE(jstart_))-ENERGY(IndexE(istart_))
      end if
      if (abs(EDIFF_) <= 1.0e-8_wp) cycle

      if ((JSTART == 1) .and. (EDIFF_ < Zero)) cycle

      IJSO = IJSO+1

      ! The energy difference is used to define the norm of the wave vector.

      rkNorm = abs(EDIFF_)/c_in_au

      ! Iterate over the quadrature points.

      OscStr(:,:) = Zero

      ! Initialize output arrays

      RAW(:,:,:) = Zero
      Aux(:,:) = Zero

      do iQuad=1,nQuad
        iVec_ = (iVec-1)*nQuad+iQuad

        ! Read or generate the wavevector

        ! Generate the wavevector associated with this quadrature
        ! point and pick up the associated quadrature weight.

        UK(:) = Rquad(1:3,iQuad)
        Wavevector(:) = rkNorm*UK(:)

        ! Note that the weights are normalized to integrate to
        ! 4*pi over the solid angles.

        Weight = Rquad(4,iQuad)
        if (.not. Do_SK) Weight = Weight/(Four*Pi)

        ! Generate the polarization vector

        if (Do_Pol) then
          pol_Vector(:,iVec_) = e_Vector-DDot_(3,UK,1,e_Vector,1)*UK
          rNorm = DDot_(3,pol_Vector(:,iVec_),1,pol_Vector(:,iVec_),1)
          if (rNorm > 1.0e-12_wp) then
            pol_Vector(:,iVec_) = pol_Vector(:,iVec_)/sqrt(rNorm)
          else
            pol_Vector(:,iVec_) = Zero
          end if
        end if

        ! Generate the property integrals associated with this
        ! direction of the wave vector k.

        iOpt = 1
        call TMOMInt(Wavevector,iOpt)

        ! Compute the transition property of the property
        ! integrals between the two states.

        ij_ = 0
        do i_=istart_,iend_
          I = IndexE(I_)
          do j_=jstart_,jend_
            J = IndexE(J_)
            EDIFF = ENERGY(J)-ENERGY(I)
            ij_ = ij_+1
            ! COMBINED SYMMETRY OF STATES:
            JOB1 = JBNUM(I)
            JOB2 = JBNUM(J)
            LSYM1 = IRREP(JOB1)
            LSYM2 = IRREP(JOB2)
            ISY12 = MUL(LSYM1,LSYM2)
            ! THE SYMMETRY CHECK MASK:
            MASK = 2**(ISY12-1)
            ! ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
            ! FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF TDMSCR
            call mk_IOFF(IOFF,nIrrep,NBASF,ISY12)
            ! CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
            ! AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES

            ! Pick up the transition density between the two states from
            ! disc. Generated in GTDMCTL.

            ISTATE = max(i,j)
            JSTATE = min(i,j)
            ij = ISTATE*(ISTATE-1)/2+JSTATE
            if (Diagonal) then
              IDISK = iDisk_TDM(I,J,1)
              iEmpty = iDisk_TDM(I,J,2)
              iOpt = 2
              iGo = 5
              call dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
              call MK_TWDM(nIrrep,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,IOFF,NBASF,ISY12)
              do IPRP=1,12
                IPROP = IPRTMOM(IPRP)
                ITYPE = 0
                if (PTYPE(IPROP) == 'HERMSING') ITYPE = 1
                if (PTYPE(IPROP) == 'ANTISING') ITYPE = 2
                if (PTYPE(IPROP) == 'HERMTRIP') ITYPE = 3
                if (PTYPE(IPROP) == 'ANTITRIP') ITYPE = 4
                LABEL = PNAME(IPROP)
                call MK_PROP(PROP,IPROP,I,J,LABEL,ITYPE,IP,NIP,SCR,nSCR,MASK,ISY12,IOFF)
              end do ! IPRP
            else

              do IPRP=1,12
                Prop(:,:,IPRTMOM(IPRP)) = Zero
              end do
              do k_=1,nState
                k = IndexE(k_)
                JOB3 = JBNUM(k)
                LSYM3 = IRREP(JOB3)
                if ((abs(EigVec(k,I)) < 1.0e-10_wp) .and. (abs(EigVec(k,J)) < 1.0e-10_wp)) cycle
                do l_=1,k
                  l = IndexE(l_)
                  JOB4 = JBNUM(l)
                  LSYM4 = IRREP(JOB4)
                  if ((abs(EigVec(l,I)) < 1.0e-10_wp) .and. (abs(EigVec(l,J)) < 1.0e-10_wp)) cycle

                  ISY34 = MUL(LSYM3,LSYM4)

                  MASK34 = 2**(ISY34-1)
                  call mk_IOFF(IOFF,nIrrep,NBASF,ISY34)

                  IDISK = iDisk_TDM(k,l,1)
                  iEmpty = iDisk_TDM(k,l,2)
                  iOpt = 2
                  iGo = 5
                  call dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,k,l)
                  call MK_TWDM(nIrrep,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,IOFF,NBASF,ISY34)

                  do IPRP=1,12
                    IPROP = IPRTMOM(IPRP)
                    ITYPE = 0
                    if (PTYPE(IPROP) == 'HERMSING') ITYPE = 1
                    if (PTYPE(IPROP) == 'ANTISING') ITYPE = 2
                    if (PTYPE(IPROP) == 'HERMTRIP') ITYPE = 3
                    if (PTYPE(IPROP) == 'ANTITRIP') ITYPE = 4
                    LABEL = PNAME(IPROP)
                    call MK_PROP(PROP,IPROP,K,L,LABEL,ITYPE,IP,NIP,SCR,nSCR,MASK34,ISY34,IOFF)
                  end do ! IPRP
                end do
              end do

              ! Transform to the new basis. Do it just for the
              ! elements we will use.

              call mma_allocate(SCR1,NSTATE,Label='SCR1')
              SCR1(:) = Zero
              do IPRP=1,12
                IPROP = IPRTMOM(IPRP)
                call DGEMM_('N','N',NSTATE,1,NSTATE,One,PROP(1,1,IPROP),NSTATE,EIGVEC(1,J),NSTATE,Zero,SCR1,NSTATE)
                PROP(I,J,IPROP) = DDot_(NSTATE,SCR1,1,EIGVEC(1,I),1)
              end do
              call mma_deallocate(SCR1)

            end if

            ! (1) the oam part

            ! The contribution to the generalized momentum operator.
            ! Note that the integrals contain nabla, but we need p,
            ! so we multiply by -i

            do iCar=1,3
              TM_R(iCar) = PROP(I,J,IPRTMOM(3+iCar))+PROP(I,J,IPRTMOM(9+iCar))
              TM_I(iCar) = -PROP(I,J,IPRTMOM(0+iCar))-PROP(I,J,IPRTMOM(6+iCar))
            end do

            ! (2) the magnetic-spin part

            ! Well the B.S term is overkill, get rid of it.
            !  Why do it when we don't do the L.S-term!

            ! Finally, evaluate the transition moment from the two
            ! different contributions.

#           ifdef _HDF5_
            ! Fix the triangular index because we are not storing the diagonal
            IJSF = IJ-ISTATE+1
            Storage(ip_w,iQuad,IJSF,iVec) = Weight
            call DCopy_(3,Wavevector,1,Storage(ip_kvector,iQuad,IJSF,iVec),1)
            call DCopy_(3,TM_R,1,Storage(ip_TMR,iQuad,IJSF,iVec),1)
            call DCopy_(3,TM_I,1,Storage(ip_TMI,iQuad,IJSF,iVec),1)
#           endif

            ! Project out the k direction from the real and imaginary components

            call DaXpY_(3,-DDot_(3,TM_R,1,UK,1),UK,1,TM_R,1)
            call DaXpY_(3,-DDot_(3,TM_I,1,UK,1),UK,1,TM_I,1)

            ! Implicitly integrate over all directions of the
            ! polarization vector to get the average value.

            TM1 = DDot_(3,TM_R,1,TM_R,1)
            TM2 = DDot_(3,TM_I,1,TM_I,1)
            TM_2 = Half*(TM1+TM2)

            ! Compute maximum and minimum oscillator strengths
            ! and the corresponding polarization vectors

            if (Do_SK) then
              TM3 = DDot_(3,TM_R,1,TM_I,1)
              Rng = sqrt((TM1-TM2)**2+Four*TM3**2)
              Aux(1,ij_) = TM_2+Half*Rng
              Aux(5,ij_) = TM_2-Half*Rng
              ! The direction for the maximum
              Ang = Half*atan2(Two*TM3,TM1-TM2)
              call daXpY_(3,cos(Ang),TM_R,1,Aux(2,ij_),1)
              call daXpY_(3,sin(Ang),TM_I,1,Aux(2,ij_),1)
              ! Normalize and compute the direction for the minimum
              ! as a cross product with k
              rNorm = DDot_(3,Aux(2,ij_),1,Aux(2,ij_),1)
              if (rNorm > 1.0e-12_wp) then
                call dScal_(3,One/sqrt(rNorm),Aux(2,ij_),1)
                Aux(6,ij_) = Aux(3,ij_)*UK(3)-Aux(4,ij_)*UK(2)
                Aux(7,ij_) = Aux(4,ij_)*UK(1)-Aux(2,ij_)*UK(3)
                Aux(8,ij_) = Aux(2,ij_)*UK(2)-Aux(3,ij_)*UK(1)
                rNorm = DDot_(3,Aux(6,ij_),1,Aux(6,ij_),1)
                call dScal_(3,One/sqrt(rNorm),Aux(6,ij_),1)
              else
                call dCopy_(3,[Zero],0,Aux(2,ij_),1)
                call dCopy_(3,[Zero],0,Aux(6,ij_),1)
              end if
            end if

            ! Oscillator strength for a specific polarization vector

            if (Do_Pol) then
              TM1 = DDot_(3,TM_R,1,pol_Vector(1,iVec_),1)
              TM2 = DDot_(3,TM_I,1,pol_Vector(1,iVec_),1)
              TM_2 = TM1*TM1+TM2*TM2
            end if

            ! Compute the oscillator strength

            F_Temp = Two*TM_2/EDIFF
            if (Do_SK) then
              Aux(1,ij_) = Two*Aux(1,ij_)/EDIFF
              Aux(5,ij_) = Two*Aux(5,ij_)/EDIFF
            end if

            ! Compute the rotatory strength

            TM_C(1) = TM_R(2)*TM_I(3)-TM_R(3)*TM_I(2)
            TM_C(2) = TM_R(3)*TM_I(1)-TM_R(1)*TM_I(3)
            TM_C(3) = TM_R(1)*TM_I(2)-TM_R(2)*TM_I(1)
            TM_2 = Two*DDot_(3,TM_C,1,UK,1)

            ! R = 3/4 * c*hbar^2/DeltaE^2 * (|T^L|^2 - |T^R|^2)

            R_Temp = 0.75_wp*c_in_au/EDIFF**2*TM_2

            ! Now let's convert this to reduced rotational strength
            ! (units of 1e-2 debye*Bohr_magneton)

            R_Temp = R_Temp*AU2REDR

            ! Save the raw oscillator and rotatory strengths in a given direction

            RAW(IQUAD,1,ij_) = F_Temp
            RAW(IQUAD,2,ij_) = R_Temp

            ! Save the direction and weight too

            RAW(IQUAD,3,ij_) = UK(1)
            RAW(IQUAD,4,ij_) = UK(2)
            RAW(IQUAD,5,ij_) = UK(3)
            RAW(IQUAD,6,ij_) = Weight

            ! Compute the oscillator and rotatory strength

            OscStr(1,ij_) = OscStr(1,ij_)+Weight*F_Temp
            OscStr(2,ij_) = OscStr(2,ij_)+Weight*R_Temp
          end do ! j_
        end do ! i_

      end do ! iQuad

      ij_ = 0
      do i_=istart_,iend_
        I = IndexE(I_)
        do j_=jstart_,jend_
          J = IndexE(J_)
          ij_ = ij_+1

          F = OscStr(1,ij_)
          R = OscStr(2,ij_)

          call Add_Info('ITMS(SF)',[F],1,6)
          call Add_Info('ROTS(SF)',[R],1,4)

          if (Do_Pol) then
            F_CHECK = abs(Aux(1,ij_))
            R_CHECK = Zero ! dummy assign
          else
            F_CHECK = abs(F)
            R_CHECK = abs(R)
          end if
          if ((F_CHECK < OSTHR) .and. (R_CHECK < RSTHR)) cycle
          A = (AFACTOR*EDIFF**2)*F

          if (iPrint == 0) then
            write(u6,*)
            if (Do_SK) then
              call CollapseOutput(1,'Transition moment strengths (spin-free states):')
              write(u6,'(3X,A)') '-----------------------------------------------'
              if (Do_Pol) then
                iVec_ = (iVec-1)*nQuad+1
                write(u6,'(4x,a,3F8.4)') 'Direction of the polarization: ',(pol_vector(k,iVec_),k=1,3)
              else
                write(u6,'(4x,a)') 'The oscillator strength is integrated over all directions of the polarization vector'
              end if
              write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3)
            else
              call CollapseOutput(1,'Isotropic transition moment strengths (spin-free states):')
              write(u6,'(3X,A)') '---------------------------------------------------------'
            end if
            if (OSTHR > Zero) write(u6,45) 'For osc. strength at least',OSTHR,'and red. rot. strength  at least',RSTHR
            write(u6,*)
            if (.not. Do_SK) then
              write(u6,'(4x,a,I4,a)') 'Integrated over ',nQuad,' directions of the wave vector'
              write(u6,'(4x,a)') 'The oscillator strength is integrated over all directions of the polarization vector'
              write(u6,*)
            end if
            write(u6,39) 'From','To','Osc. strength','Red. rot. str.','Total A (sec-1)'
            write(u6,40)
            iPrint = 1
          end if

          ! Regular print

          if (F_CHECK < OSTHR) then
            write(u6,46) I,J,'below threshold',R,A
          else if (R_CHECK < RSTHR) then
            ! Don't print rot. str. if below threshold
            write(u6,47) I,J,F,'below threshold',A
          else
            write(u6,33) I,J,F,R,A
          end if

          if (Do_SK) then
            write(u6,50) 'maximum',Aux(1,ij_),'for polarization direction:',Aux(2,ij_),Aux(3,ij_),Aux(4,ij_)
            write(u6,50) 'minimum',Aux(5,ij_),'for polarization direction:',Aux(6,ij_),Aux(7,ij_),Aux(8,ij_)
          end if

          ! Printing raw (unweighted) and direction for every transition

          if (PRRAW) then
            write(u6,*)
            write(u6,*)
            write(u6,41) 'From','To','Raw osc. str.','Rot. str.','kx','ky','kz'
            write(u6,32)
            do IQUAD=1,NQUAD
              write(u6,33) I,J,RAW(IQUAD,1,ij_),RAW(IQUAD,2,ij_),RAW(IQUAD,3,ij_),RAW(IQUAD,4,ij_),RAW(IQUAD,5,ij_)
            end do
            write(u6,32)
            write(u6,*)
          end if

          ! Printing weighted and direction for every transition

          if (PRWEIGHT) then
            write(u6,*)
            write(u6,*)
            write(u6,41) 'From','To','Weig. osc. str.','Rot. str.','kx','ky','kz'
            write(u6,32)
            do IQUAD=1,NQUAD
              Weight = RAW(IQUAD,6,ij_)
              write(u6,33) I,J,RAW(IQUAD,1,ij_)*Weight,RAW(IQUAD,2,ij_)*Weight,RAW(IQUAD,3,ij_),RAW(IQUAD,4,ij_),RAW(IQUAD,5,ij_)
            end do
            write(u6,32)
            write(u6,*)
          end if
        end do
      end do

    end do
  end do

  if (iPrint == 1) then
    write(u6,40)
    if (Do_SK) then
      call CollapseOutput(0,'Transition moment strengths (spin-free states):')
    else
      call CollapseOutput(0,'Isotropic transition moment strengths (spin-free states):')
    end if
  end if

end do ! iVec

#ifdef _HDF5_
call mh5_put_dset(wfn_sfs_tm,pack(Storage,.true.))
call mma_deallocate(Storage)
#endif

! Deallocate some arrays.

call mma_deallocate(RAW)
call mma_deallocate(IP)
call mma_deallocate(OscStr)
call mma_deallocate(Aux)
if (TMOgroup) then
  call mma_DeAllocate(TMOgrp1)
  call mma_DeAllocate(TMOgrp2)
end if
if (Do_Pol) call mma_deallocate(pol_Vector)
call mma_deallocate(SCR)
call mma_deAllocate(TDMZZ)
call mma_deAllocate(TSDMZZ)
call mma_deAllocate(WDMZZ)

#ifdef _TIME_TMOM_
call CWTime(TCpu2,TWall2)
write(u6,*) 'Time for TMOM (SF) ',TCpu2-TCpu1,TWall2-TWall1
#endif

! Do some cleanup

call mma_deAllocate(Rquad)
call ClsSew()
!***********************************************************************
!                                                                      *
!     End of section for transition moments                            *
!                                                                      *
!***********************************************************************

call FinishUp()

222 format(5X,2(1X,I4),5X,3(1X,ES18.8))
30 format(5X,A,1X,ES15.8)
31 format(5X,2(1X,A4),6X,A15,1X,A47,1X,A15)
32 format(5X,95('-'))
33 format(5X,2(1X,I4),5X,5(1X,ES15.8))
34 format(5X,2(1X,A4),5X,4(1X,A15),1X,A)
35 format(5X,31('-'))
36 format(5X,2(1X,I4),6X,15('-'),1X,ES15.8,1X,A15)
37 format(5X,2(1X,I4),6X,15('-'),1X,A15,1X,ES15.8)
38 format(5X,2(1X,I4),6X,F15.6,4(1X,ES15.8))
39 format(5X,2(1X,A4),5X,3(1X,A15))
40 format(5X,63('-'))
41 format(5X,2(1X,A4),5X,5(1X,A15))
42 format(5X,79('-'))
43 format(12X,A8,6(1X,ES15.8))
44 format(20X,6(1X,A15))
45 format(4X,2(A,1X,ES15.8,1X))
46 format(5X,2(1X,I4),5X,(1X,A15),2(1X,ES15.8))
47 format(5X,2(1X,I4),5X,(1X,ES15.8),(1X,A15),(1X,ES15.8))
49 format(5X,A,1X,ES15.8,1X,A)
50 format(10X,A7,3X,1(1X,ES15.8),5X,A27,3(1X,F7.4))

contains

subroutine FinishUp()

# ifdef _DEBUGPRINT_
  write(u6,*) 'end of eigctl: BLUBB debug print of property matrix'
  do istate=1,nstate
    do jstate=1,nstate
      do IPROP=1,NPROP
        if (abs(prop(istate,jstate,iprop)) > 1.0e-14_wp) &
          write(u6,*) 'prop(',istate,',',jstate,',',iprop,') = ',prop(istate,jstate,iprop)
      end do
    end do
  end do
# endif
  call mma_DeAllocate(IndexE)

end subroutine FinishUp

end subroutine EigCtl
