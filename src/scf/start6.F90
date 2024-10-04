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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************

subroutine Start6(FName,LuOrb,CMO,mBB,nD,EOrb,OccNo,mmB)
!***********************************************************************
!                                                                      *
!     purpose: Generate constrained orbitals from INPORB               *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use SpinAV, only: Do_SpinAV, DSC
use InfSCF, only: DoCholesky, E_nondyn, Erest_xc, FileOrb_id, IndxC, isHDF5, MaxBas, MxConstr, nBas, nBB, nBT, nConstr, nDel, &
                  nFro, nnB, nOcc, nOrb, nSym, s2CNO, VTitle
use Cholesky, only: ChFracMem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: FName
integer(kind=iwp), intent(in) :: LuOrb, mBB, nD, mmB
real(kind=wp), intent(out) :: CMO(mBB,nD), EOrb(mmB,nD), OccNo(mmB,nD)
integer(kind=iwp) :: i, ibas, ic1, ic2, iCMO, iComp, iDaa, iDbb, iDSc, iDummy(1), iErr, Indx, iOcc, iOff, iOpt, ipDaa, ipDbb, &
                     ipDScc, ipMK, ipML, iRC, iSym, iSymLbl, iWFType, j, jc, ji, jOcc, jOff, k, kc, kc1, kc2, kDSc, kk, kkc, kks, &
                     l, lc, lc1, lc2, llc, lls, lOcc, lOff, lsq, ltri, Lu_, mAdCMOO, mOff, nBD(8), nDiff_ab, nHoles(8), nIF(8), &
                     nRASO(8), nSsh(8), nSsh_ab(8), nTmp(8), nZero(8)
real(kind=wp) :: Dummy(1), ThrD, xNorm, xOkk, yOkk
character(len=62) :: Line
character(len=8) :: Label
integer(kind=iwp), allocatable :: ID_vir(:), IndT(:), Match(:,:)
real(kind=wp), allocatable :: Corb(:), Da(:,:), SAV(:), SLT(:), SQ(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

Erest_xc = Zero

if (.not. DoCholesky) then
  write(u6,*)
  write(u6,*) ' ERROR in Constrained SCF: problem in start6.'
  write(u6,*) '*** Constrained NOs implemented only with CD or RI.'
  write(u6,*) '*** Use Cholesky or RICD in Seward and rerun! *****'
  call Abend()
end if

do iSym=1,nSym
  nDiff_ab = nOcc(iSym,1)-nOcc(iSym,2)
  if (nDiff_ab < 0) then
    write(u6,*)
    write(u6,*) ' ERROR in Constrained SCF: problem in start6.'
    write(u6,*) '*** #alpha < #beta not permitted in CNOs    ***'
    write(u6,*) '*** Change SCF input accordingly and rerun! ***'
    call Abend()
  end if
  nHoles(iSym) = nDiff_ab
end do

write(u6,*) ' -------------------------------------------------'
if (Do_SpinAV) then
  write(u6,*) ' Spin-averaged wavelets (+/-) '
else
  write(u6,*) ' Configuration of the constrained spins (up/down) '
end if
write(u6,*) ' -------------------------------------------------'
do iSym=1,nSym
  write(u6,'(1X,A,I1)') ' sym: ',iSym
  Line(1:14) = '         (+) '
  k = 15
  do j=1,nConstr(iSym)
    if (indxC(j,1,iSym) == 1) then
      if (Do_SpinAV) then
        Line(k:k+2) = ' + '
      else
        Line(k:k+2) = ' u '
      end if
    else if (indxC(j,1,iSym) == 2) then
      if (Do_SpinAV) then
        Line(k:k+2) = ' - '
      else
        Line(k:k+2) = ' d '
      end if
    else
      Line(k:k+2) = '   '
    end if
    k = k+3
  end do
  write(u6,*) Line(1:k-1)
  Line(1:14) = '         (-) '
  k = 15
  do j=1,nConstr(iSym)
    if (indxC(j,2,iSym) == 1) then
      if (Do_SpinAV) then
        Line(k:k+2) = ' + '
      else
        Line(k:k+2) = ' u '
      end if
    else if (indxC(j,2,iSym) == 2) then
      if (Do_SpinAV) then
        Line(k:k+2) = ' - '
      else
        Line(k:k+2) = ' d '
      end if
    else
      Line(k:k+2) = '   '
    end if
    k = k+3
  end do
  write(u6,*) Line(1:k-1)
end do
write(u6,*) ' -------------------------------------------------'

Lu_ = LuOrb
call mma_Allocate(IndT,nnB,Label='IndT')
if (isHDF5) then
  call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,CMO,OccNo,EOrb,IndT)
else
  call RdVec_(FName,Lu_,'COEI',0,nSym,nBas,nOrb,CMO,Dummy,OccNo,Dummy,EOrb,Dummy,IndT,VTitle,1,iErr,iWFtype)
end if
call RdTwoEnrg(Lu_,E_nondyn)
call VecSort(nSym,nBas,nBas,CMO,OccNo,IndT,0,iDummy,iErr)
indx = 1
nZero(1:nSym) = 0
nTmp(1:nSym) = 0
nIF(1:nSym) = 0
nRASO(1:nSym) = 0
do iSym=1,nSym
  nDiff_ab = 0
  do iBas=1,nBas(iSym)
    select case (IndT(indx))
      case (1,2)
        nIF(iSym) = nIF(iSym)+1 ! froz + inac orbitals
      case (3)
        nIF(iSym) = nIF(iSym)+1  ! electrons (place them in RAS1)
        nDiff_ab = nDiff_ab+1
      case (4,5)
        nRASO(iSym) = nRASO(iSym)+1
      case (7)
        nTmp(iSym) = nTmp(iSym)+1
    end select
    indx = indx+1
  end do
  if (nRASO(iSym) /= 2*nConstr(iSym)) then
    write(u6,*) ' ERROR in Constrained SCF: problem in start6.'
    write(u6,*) ' Detected inconsistency between # of partially occupied orbitals and # of constraints. Sym: ',iSym
    call Abend()
  end if
  if (nHoles(iSym) /= nDiff_ab) then
    write(u6,*) ' ERROR in Constrained SCF: problem in start6.'
    write(u6,*) ' Detected inconsistency between # of excess alpha orbitals and # of RAS1 orbitals. Sym: ',iSym
    call Abend()
  end if
  if (nOrb(iSym) > nBas(iSym)-nTmp(iSym)) then
    nOrb(iSym) = nBas(iSym)-nTmp(iSym)
    nDel(iSym) = nTmp(iSym)
  end if
end do

call TrimCMO(CMO,nSym,nBas,nOrb)
call TrimEor(EOrb,nSym,nBas,nOrb)
call mma_deallocate(IndT)

call Setup_SCF()

nBD(1) = 0
do iSym=2,nSym
  nBD(iSym) = nBD(iSym-1)+nTri_Elem(nBas(iSym-1))
end do
call mma_allocate(Da,nBT,2,Label='Da')
Da(:,:) = Zero
call mma_allocate(Match,2,MxConstr,Label='Match')
call mma_allocate(Corb,MaxBas,Label='Corb')

if (Do_SpinAV) call mma_allocate(SAV,2*MaxBas**2,Label='SAV')

iOff = 1
jOff = 0
do iSym=1,nSym
  CMO(iOff:iOff+nBas(iSym)*nOrb(iSym)-1,2) = CMO(iOff:iOff+nBas(iSym)*nOrb(iSym)-1,1)
  EOrb(jOff+1:jOff+nOrb(iSym),2) = EOrb(jOff+1:jOff+nOrb(iSym),1)
  lOcc = 1+jOff+nIF(iSym)
  OccNo(1:nRASO(iSym),2) = OccNo(lOcc:lOcc+nRASO(iSym)-1,1)
  call BestMatch(nConstr(iSym),nRASO(iSym),OccNo(1,2),Match,MxConstr)
  do i=1,nConstr(iSym)
    k = Match(1,i) ! (+) wavelet
    jOcc = jOff+nIF(iSym)+k
    xOkk = OccNo(jOcc,1)/Two
    kc = iOff+nBas(iSym)*(nIF(iSym)+k-1)
    l = Match(2,i)  ! (-) wavelet
    iOcc = jOff+nIF(iSym)+l
    yOkk = OccNo(iOcc,1)/Two
    xnorm = sqrt(abs(xOkk)+abs(yOkk)) ! ensures correct normaliz
    lc = iOff+nBas(iSym)*(nIF(iSym)+l-1)
    xOkk = sqrt(abs(xOkk))/xnorm
    yOkk = sqrt(abs(yOkk))/xnorm
    if (Do_SpinAV) then
      kkc = 1+nBas(iSym)*(k-1)
      llc = 1+nBas(iSym)*(nConstr(iSym)+l-1)
      SAV(kkc:kkc+nBas(iSym)-1) = CMO(kc:kc+nBas(iSym)-1,1)-xOkk*CMO(lc:lc+nBas(iSym)-1,1)
      SAV(llc:llc+nBas(iSym)-1) = xOkk*CMO(lc:lc+nBas(iSym)-1,1)+yOkk*CMO(kc:kc+nBas(iSym)-1,1)
    end if
    Corb(1:nBas(iSym)) = xOkk*CMO(kc:kc+nBas(iSym)-1,1)
    CMO(kc:kc+nBas(iSym)-1,1) = Corb(1:nBas(iSym))+yOkk*CMO(lc:lc+nBas(iSym)-1,1)
    CMO(lc:lc+nBas(iSym)-1,1) = Corb(1:nBas(iSym))-yOkk*CMO(lc:lc+nBas(iSym)-1,1)
  end do
  jc = 1
  kc = nConstr(iSym)+1
  do i=1,nConstr(iSym)
    l = Match(indxC(i,2,iSym),i)
    lc1 = iOff+nBas(iSym)*(nIF(iSym)+l-1)
    lc2 = iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+jc-1)
    CMO(lc2:lc2+nBas(iSym)-1,2) = CMO(lc1:lc1+nBas(iSym)-1,1)
    k = Match(indxC(i,1,iSym),i)
    kc1 = iOff+nBas(iSym)*(nIF(iSym)+k-1)
    kc2 = iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+kc-1)
    CMO(kc2:kc2+nBas(iSym)-1,2) = CMO(kc1:kc1+nBas(iSym)-1,1)
    jc = jc+1
    kc = kc+1
  end do
  kc = nConstr(iSym)+1
  do i=1,nConstr(iSym)
    ic1 = iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+i-1)
    ic2 = iOff+nBas(iSym)*(nIF(iSym)+kc-1)
    CMO(ic2:ic2+nBas(iSym)-1,1) = CMO(ic1:ic1+nBas(iSym)-1,2)
    kc1 = iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+kc-1)
    kc2 = iOff+nBas(iSym)*(nIF(iSym)+i-1)
    CMO(kc2:kc2+nBas(iSym)-1,1) = CMO(kc1:kc1+nBas(iSym)-1,2)
    kc = kc+1
  end do
  kc = nConstr(iSym)+1  ! wavelets in virt space
  if (Do_SpinAV) then
    do i=1,nConstr(iSym)
      k = Match(1,i)
      kks = 1+nBas(iSym)*(k-1)
      l = Match(2,i)
      lls = 1+nBas(iSym)*(nConstr(iSym)+l-1)
      mOff = iOff+nBas(iSym)*(nIF(iSym)+kc-1)
      kk = indxC(i,1,iSym)
      if (kk == 1) then ! => (+) wavelet is in alpha
        ipMK = mOff
        ipML = mOff-nBas(iSym)*nHoles(iSym)
        CMO(ipMK:ipMK+nBas(iSym)-1,1) = SAV(kks:kks+nBas(iSym)-1)
        CMO(ipML:ipML+nBas(iSym)-1,2) = SAV(lls:lls+nBas(iSym)-1)
      else if (kk == 2) then
        ipMK = mOff-nBas(iSym)*nHoles(iSym)
        ipML = mOff
        CMO(ipMK:ipMK+nBas(iSym)-1,2) = SAV(kks:kks+nBas(iSym)-1)
        CMO(ipML:ipML+nBas(iSym)-1,1) = SAV(lls:lls+nBas(iSym)-1)
      else
        ipMK = 666666  ! avoid compiler wrngs
        ipML = 666666
        write(u6,*) ' Start6: wrong indxC value: ',kk
        call Abend()
      end if
      kc = kc+1
    end do
  end if
  iOff = iOff+nBas(iSym)*nOrb(iSym)
  jOff = jOff+nOrb(iSym)
end do

if (Do_SpinAV) then
  call mma_deallocate(SAV)
  call mma_Allocate(DSc,nBB,Label='DSc')
  DSC(:) = Zero
end if

iOff = 1
lOff = 0
do iSym=1,nSym
  ipDaa = 1+nBD(iSym)
  mAdCMOO = iOff+nBas(iSym)*nIF(iSym)
  call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym), &
                 One,CMO(mAdCMOO,1),nBas(iSym), &
                 CMO(mAdCMOO,1),nBas(iSym), &
                 Zero,Da(ipDaa,1),nBas(iSym))
  ipDbb = 1+nBD(iSym)
  mAdCMOO = iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym))
  call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym), &
                 One,CMO(mAdCMOO,2),nBas(iSym), &
                 CMO(mAdCMOO,2),nBas(iSym), &
                 Zero,Da(ipDbb,2),nBas(iSym))

  if (Do_SpinAV) then
    ipDScc = lOff
    do j=1,nBas(iSym)
      do i=1,j
        ji = iTri(j,i)
        iDaa = ipDaa-1+ji
        iDbb = ipDbb-1+ji
        iDSc = ipDScc+nBas(iSym)*(j-1)+i
        DSc(iDSc) = Half*(Da(iDaa,1)-Da(iDbb,2))
        kDSc = ipDScc+nBas(iSym)*(i-1)+j
        DSc(kDSc) = DSc(iDSc)
      end do
    end do
    lOff = lOff+nBas(iSym)**2
  end if

  do j=1,nBas(iSym)
    do i=1,j-1
      ji = iTri(j,i)
      iDaa = ipDaa-1+ji
      Da(iDaa,1) = Two*Da(iDaa,1)
      iDbb = ipDbb-1+ji
      Da(iDbb,2) = Two*Da(iDbb,2)
    end do
  end do
  iOff = iOff+nBas(iSym)*nOrb(iSym)
end do

call Cho_X_init(irc,ChFracMem)
if (irc /= 0) then
  call WarningMessage(2,'Start6. Non-zero rc in Cho_X_init.')
  call Abend()
end if
!----------------------------------------------------------------------*
call Get_Fmat_nondyn(Da(:,1),Da(:,2),nBT,.false.)
!----------------------------------------------------------------------*

call Cho_X_Final(irc)
if (irc /= 0) then
  call WarningMessage(2,'Start6. Non-zero rc in Cho_X_Final.')
  call Abend()
end if

call mma_deallocate(Da)
call mma_deallocate(Corb)
call mma_deallocate(Match)

iOff = 0
jOff = 0
do iSym=1,nSym
  OccNo(iOff+1:iOff+nOcc(iSym,1),1) = One
  OccNo(iOff+nOcc(iSym,1)+1:iOff+nOrb(iSym),1) = Zero
  OccNo(iOff+1:iOff+nOcc(iSym,2),2) = One
  OccNo(iOff+nOcc(iSym,2)+1:iOff+nOrb(iSym),2) = Zero
  iOff = iOff+nOrb(iSym)
end do

call mma_allocate(SLT,nBT,Label='SLT')
isymlbl = 1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'Mltpl  0'
iComp = 1
call RdOne(irc,iOpt,Label,iComp,SLT,isymlbl)
if (irc /= 0) then
  write(u6,*) ' Start6 : error in getting overlap matrix '
  call Abend()
end if
call s2calc(CMO(:,1),CMO(:,2),SLT,nOcc(:,1),nOcc(:,2),nBas,nOrb,nSym,s2CNO)

if (.not. Do_SpinAV) then
  write(u6,'(A,f9.6)') '  Initial value of Total Spin, S(S+1): ',s2CNO
  write(u6,*) ' -------------------------------------------------'
end if
write(u6,*)

!----------------------------------------------------------------------*
!  Virtual space must be orthogonal to the occupied space              *
!----------------------------------------------------------------------*
if (Do_SpinAV) then
  nOcc(1:nSym,1) = nOcc(1:nSym,1)+nConstr(1:nSym)
  nOcc(1:nSym,2) = nOcc(1:nSym,2)+nConstr(1:nSym)
end if

Thrd = 1.0e-6_wp
nSsh(1:nSym) = nOrb(1:nSym)-nOcc(1:nSym,1)-nFro(1:nSym)
nSsh_ab(1:nSym) = nOrb(1:nSym)-nOcc(1:nSym,2)-nFro(1:nSym)

call mma_allocate(SQ,nBB,Label='SQ')
ltri = 1
lsq = 1
do iSym=1,nSym
  call Square(SLT(ltri),SQ(lsq),1,nBas(iSym),nBas(iSym))
  ltri = ltri+nTri_Elem(nBas(iSym))
  lsq = lsq+nBas(iSym)**2
end do
call mma_allocate(ID_vir,nnB,Label='ID_vir')
call Cho_ov_Loc(irc,Thrd,nSym,nBas,nOcc(:,1),nZero,nZero,nSsh,CMO(:,1),SQ,ID_vir)

if (irc /= 0) then
  write(u6,*) ' Start6 : error in getting alpha virt MOs '
  call Abend()
end if

iOff = 1
do iSym=1,nSym
  iCMO = iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym,1))
  call Ortho_Orb(CMO(iCMO:,1),SQ,nBas(iSym),nSsh(iSym),2,.true.)
  iOff = iOff+nBas(iSym)*nOrb(iSym)
end do
!call ChkOrt(1,Whatever) ! silent

call Cho_ov_Loc(irc,Thrd,nSym,nBas,nOcc(1,2),nZero,nZero,nSsh_ab,CMO(:,2),SQ,iD_vir)

if (irc /= 0) then
  write(u6,*) ' Start6 : error in getting beta virt MOs '
  call Abend()
end if

iOff = 1
do iSym=1,nSym
  iCMO = iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym,2))
  call Ortho_Orb(CMO(iCMO:,2),SQ,nBas(iSym),nSsh_ab(iSym),2,.true.)
  iOff = iOff+nBas(iSym)*nOrb(iSym)
end do
!call ChkOrt(2,Whatever) ! silent

if (Do_SpinAV) then ! reset # of occupied
  nOcc(1:nSym,1) = nOcc(1:nSym,1)-nConstr(1:nSym)
  nOcc(1:nSym,2) = nOcc(1:nSym,2)-nConstr(1:nSym)
end if

call mma_deallocate(ID_vir)
call mma_deallocate(SQ)
call mma_deallocate(SLT)

return

end subroutine Start6
