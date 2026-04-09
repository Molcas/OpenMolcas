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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               2016,2017, Roland Lindh                                *
!***********************************************************************

subroutine FinalSCF()
!***********************************************************************
!                                                                      *
!     purpose: perform final calculations                              *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Interfaces_SCF, only: dOne_SCF
use OFembed, only: Do_OFemb, FMaux, NDSD
use SpinAV, only: DSc
use InfSCF, only: BName, CMO, Darwin, Dens, DMOMax, Do_Addc, Do_Tw, DSCF, E1V, E2V, EneV, EOrb, Falcon, FMOMax, FockAO, FThr, &
                  iPrForm, iPrint, iStatPrn, kIVO, KntE, KSDFT, MaxBas, MaxBXO, MssVlc, NamFld, nBas, nBB, nBO, nBT, nD, nDel, &
                  nDens, nFld, nFro, nIter, nIterP, nnB, nnO, nOcc, NoProp, nOrb, nSym, OccNo, OneHam, Ovrlp, TimFld, TwoHam
use SCFFiles, only: LuOut
#ifdef _EFP_
use EFP_Module, only: EFP_Instance
use EFP, only: EFP_Shutdown
#endif
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use SCFWFn, only: wfn_mocoef, wfn_mocoef_a, wfn_mocoef_b, wfn_occnum, wfn_occnum_a, wfn_occnum_b, wfn_orbene, wfn_orbene_a, &
                  wfn_orbene_b, wfn_tpidx, wfn_tpidx_a, wfn_tpidx_b
#endif
#ifdef _FDE_
use Embedding_Global, only: embPot, embWriteEsp
use InfSCF, only: nAtoms
#endif
use Molcas, only: MxSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iBas, iCMO, iD, iFld, iFock, iiOrb, ij, IndType(7,8), iOpt, iOrb, iRC, iRef, iSym, iSymLb, iVirt, iWFType, &
                     jBas, jRef, jVirt, kBas, kl, lk, nFldP, nOccMax(nSym), nOccMin(nSym)
real(kind=wp) :: Dummy(1), TCPU1, TCPU2, TotCpu, TWall1, TWall2
logical(kind=iwp) :: FstItr
character(len=128) :: OrbName
character(len=80) :: Note
character(len=60) :: Frmt
character(len=8) :: Method, RlxLbl, What
real(kind=wp), allocatable :: CMOn(:), DMat(:,:), E_Or(:,:), Epsn(:), Etan(:), GVFck(:,:), Scrt1(:,:), Scrt2(:,:), Temp(:)
logical(kind=iwp), external :: RF_On, Langevin_On, PCM_On
#ifdef _HDF5_
integer(kind=iwp) :: IndTypeT(8,7), nSSh(mxSym), nZero(mxSym)
character, allocatable :: typestring(:)
#endif
#ifdef _EFP_
logical(kind=iwp), external :: EFP_On
#endif

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Read remaining one-electron integrals
call R1IntB()

call CWTime(TCpu1,TWall1)

What = 'COEI'

call SorbCMOs(CMO,nBB,nD,EOrb,OccNo,nnB,nBas,nOrb,nSym)

call Put_darray('SCF orbitals',CMO(1,1),nBB)
if (nD == 2) call Put_darray('SCF orbitals_ab',CMO(1,2),nBB)

if (nIter(nIterP) <= 0) then

  FstItr = .true.
  call SCF_Energy(FstItr,E1V,E2V,EneV)
  Dens(:,:,nDens) = Dens(:,:,1)
  TwoHam(:,:,nDens) = TwoHam(:,:,1)

  DMOMax = Zero
  FMOMax = Zero

else

  ! Compute improved virtuals (if needed)
  if (kIvo == 1) then
    do iD=1,nD
      call IvoGen(OneHam,nBT,CMO(:,iD),nBO,EOrb(:,iD),nnO,nOcc(:,iD))
    end do
  end if

  !------------------------------------------------
  !----- Write stuff for gradient calculation -----
  !------------------------------------------------

  ! Add Fock operator matrix

  ! Note that we need an array with four additional
  ! elements due to the use of WrOne, which adds
  ! four elements for some auxiliary information!

  call mma_allocate(Temp,nBT+4,Label='Temp')
  Temp(:) = Zero

  do iD=1,nD
    Temp(1:nBT) = FockAO(:,iD)
    iRc = -1
    iOpt = 0
    RlxLbl = 'Fock Op '
    if (iD == 2) RlxLbl = 'Fock Op2'
    iSymLb = 1
    call WrOne(iRc,iOpt,RlxLbl,1,Temp,iSymLb)
    if (iRc /= 0) then
      write(u6,*) 'Final: Error writing on ONEINT'
      write(u6,'(A,A)') 'RlxLbl=',RlxLbl
      call Abend()
    end if
  end do
  call mma_deallocate(Temp)

  ! Construct the generalized variational Fock matrix

  call mma_allocate(GVFck,nBT,nD,Label='GVFck')
  call mma_allocate(Scrt1,MaxBas**2,nD,Label='Scrt1')
  call mma_allocate(Scrt2,MaxBxO,nD,Label='Scrt2')

  iFock = 1
  iCMo = 1
  do iSym=1,nSym
    if (nOrb(iSym) <= 0) cycle

    do iD=1,nD
      call Square(FockAO(iFock,iD),Scrt1(1,iD),1,nBas(iSym),nBas(iSym))
      ! Transform to MO basis
      call DGEMM_('T','N',nOrb(iSym),nBas(iSym),nBas(iSym), &
                  One,CMO(iCMo,iD),nBas(iSym), &
                  Scrt1(1,iD),nBas(iSym), &
                  Zero,Scrt2(1,iD),nOrb(iSym))
      call DGEMM_('N','N',nOrb(iSym),nOrb(iSym),nBas(iSym), &
                  One,Scrt2(1,iD),nOrb(iSym), &
                  CMO(iCMo,iD),nBas(iSym), &
                  Zero,Scrt1(1,iD),nOrb(iSym))

      ! Set elements with both indices virtual to zero
      do iVirt=1,nOrb(iSym)-nOcc(iSym,iD)
        jVirt = 1+nOcc(iSym,iD)*nOrb(iSym)+nOcc(iSym,iD)+(iVirt-1)*nOrb(iSym)
        Scrt1(jVirt:jVirt+nOrb(iSym)-nOcc(iSym,iD)-1,iD) = Zero
      end do
      ! Now project back to the SO basis
      call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nOrb(iSym), &
                  One,CMO(iCMo,iD),nBas(iSym), &
                  Scrt1(1,iD),nOrb(iSym), &
                  Zero,Scrt2(1,iD),nBas(iSym))
      call DGEMM_('N','T',nBas(iSym),nBas(iSym),nOrb(iSym), &
                  One,Scrt2(1,iD),nBas(iSym), &
                  CMO(iCMo,iD),nBas(iSym), &
                  Zero,Scrt1(1,iD),nBas(iSym))

      ij = iFock
      do iBas=1,nBas(iSym)
        do jBas=1,iBas-1
          kl = nBas(iSym)*(jBas-1)+iBas
          lk = nBas(iSym)*(iBas-1)+jBas
          GVFck(ij,iD) = Scrt1(kl,iD)+Scrt1(lk,iD)
          ij = ij+1
        end do
        kl = nBas(iSym)*(iBas-1)+iBas
        GVFck(ij,iD) = Scrt1(kl,iD)
        ij = ij+1
      end do

    end do ! iD

    iFock = iFock+nTri_Elem(nBas(iSym))
    iCMo = iCMo+nBas(iSym)*nOrb(iSym)

  end do  ! iSym
  call mma_deallocate(Scrt2)
  call mma_deallocate(Scrt1)

  ! Add elementary info
  Method = 'RHF-SCF '
  if (nD == 2) Method = 'UHF-SCF '

  if (kIvo /= 0) Method = 'IVO-SCF '
  if (KSDFT /= 'SCF') Method = 'KS-DFT  '
  call Put_cArray('Relax Method',Method,8)
  !call Put_Energy(EneV)
  call Store_Energies(1,[EneV],1)
  call Put_dScalar('SCF energy',EneV)
  !if (nD == 2) call Put_dScalar('Ener_ab',EneV_ab)
  call Put_iArray('nIsh',nOcc(1,1),nSym)
  if (nD == 2) call Put_iArray('nIsh_ab',nOcc(1,2),nSym)
  call Put_iArray('nOrb',nOrb,nSym)
  call Put_iArray('nDel',nDel,nSym)
  call Put_iArray('nFroPT',nFro,nSym)    ! for Cholesky-CC
  call Put_iArray('nDelPT',nDel,nSym)    !

  ! Add MO-coefficients
  call Put_dArray('Last orbitals',CMO(1,1),nBB)
  if (nD == 2) call Put_dArray('CMO_ab',CMO(1,2),nBB)

  ! Add one body density matrix in AO/SO basis

  ! If the density matrix is explicitly calculated in the beginning we will
  ! use that one instead.
  if (nIter(nIterP) <= 0) then
    call DensAB(nBT,nDens,nD,Dens)
  else
    call mma_allocate(DMat,nBT,nD,Label='DMat')

    if (nD == 1) then
      call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(:,1),nBB,OccNo(:,1),DMat(:,1),.true.)
    else
      call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(:,1),nBB,OccNo(:,1),DMat(:,1),.true.)
      call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(:,2),nBB,OccNo(:,2),DMat(:,2),.false.)
      call Put_dArray('D1aoVar',Dummy,0) ! Undefined the field.
    end if

    call DensAB(nBT,1,nD,DMat)
    call mma_deallocate(DMat)
  end if

  if (nD == 2) then
    call mma_allocate(CMOn,nBB,Label='CMOn')
    call mma_allocate(Etan,nnB,Label='Etan')
    call mma_allocate(Epsn,nnB,Label='Epsn')
    call NatoUHF(Dens(:,1,1),Dens(:,2,1),FockAO(:,1),FockAO(:,2),nBT,CMO(:,1),nBB,Ovrlp,CMOn,Etan,Epsn,nnB,nSym,nBas,nOrb)
    call PadCMO(CMOn,nSym,nBas,nOrb)
    call PadEor(Etan,nSym,nBas,nOrb)
    call PadEor(Epsn,nSym,nBas,nOrb)
  end if

  ! Add generalized Fock matrix
  if (nD == 2) then
    GVFck(:,1) = GVFck(:,1)+GVFck(:,2)
  else
    GVFck(:,1) = Two*GVFck(:,1)
  end if
  call Put_dArray('FockOcc',GVFck(:,1),nBT)
  call mma_deallocate(GVFck)

  ! Add SCF orbital energies
  call Put_dArray('OrbE',EOrb(:,1),nnB)
  if (nD == 2) call Put_dArray('OrbE_ab',EOrb(:,2),nnB)

  call Put_dScalar('Thrs',Fthr)
end if

#ifdef _FDE_
! Embedding
if (embPot .and. embWriteEsp) call embPotOutput(nAtoms,Dens)
#endif

! t.t.;
! store Fock matrix in Runfile in case of fragment calculation
if (Falcon) then
  !write(u6,*) 'Fock matrix is written in RunFile.'
  !write(u6,*) 'fck:'
  !write(u6,*) 'nbt=',nbt
  !write(u6,*) 'ndens=',ndens
  !write(u6,*) (FockAO(itt),itt=1,nbt)
  call Put_dArray('Fragment_Fock',FockAO(1,1),nBT)
end if

! t.t.; end
! Final printout
! original PrFin was splitted to 3 parts to run in UHF mode

nDel(1:nSym) = nBas(1:nSym)-nOrb(1:nSym)
call PrFin0(Dens(:,1,nDens),Dens(:,nD,nDens),nBT,EOrb(:,1),nnB,CMO(:,1),nBO,KntE)

do iD=1,nD
  call PrFin(OneHam,Ovrlp,Dens(:,iD,nDens),TwoHam(:,iD,nDens),nBT,EOrb(:,iD),OccNo(:,iD),nnB,CMO(:,iD),nBO,Note,iD-1,MssVlc,Darwin)

  call PrFin2(Ovrlp,nBT,OccNo(:,iD),nnB,CMO(:,iD),nBO,Note)
end do

if ((iPrint >= 2) .and. (nD == 2)) &
  call PriMO('Natural orbitals',.true.,.true.,Zero,2.0e9_wp,nSym,nBas,nOrb,BName,Epsn,Etan,CMOn,iPrForm)
! Calculate expectation values
if (.not. NoProp) then
  if (iPrint >= 3) write(u6,'(/6X,A)') 'Expectation values of various operators'
  call Prpt()
end if

! make a fix for energies for deleted orbitals
call mma_allocate(E_Or,nnB,nD,Label='E_Or')

iRef = 1
jRef = 1
do iSym=1,nSym
  iiOrb = nOrb(iSym)-nDel(iSym)
  do iOrb=1,iiOrb
    E_Or(iRef,1:nD) = EOrb(jRef,1:nD)
    iRef = iRef+1
    jRef = jRef+1
  end do
  do iOrb=1,nDel(iSym)
    E_Or(iRef,1:nD) = 1000.0_wp
    iRef = iRef+1
  end do
end do

if (nD == 1) then
  IndType(1,1:nSym) = nFro(1:nSym)
  IndType(2,1:nSym) = nOcc(1:nSym,1)
  IndType(3,1:nSym) = 0
  IndType(4,1:nSym) = 0
  IndType(5,1:nSym) = 0
  IndType(6,1:nSym) = nOrb(1:nSym)-nFro(1:nSym)-nOcc(1:nSym,1)-nDel(1:nSym)
  IndType(7,1:nSym) = nDel(1:nSym)
else
  nOccMax(:) = max(nOcc(1:nSym,1),nOcc(1:nSym,2))
  nOccMin(:) = min(nOcc(1:nSym,1),nOcc(1:nSym,2))
  IndType(1,1:nSym) = nFro(1:nSym)
  IndType(2,1:nSym) = nOccMin(:)
  IndType(3,1:nSym) = 0
  IndType(4,1:nSym) = nOccMax(:)-nOccMin(:)
  IndType(5,1:nSym) = 0
  IndType(6,1:nSym) = nOrb(1:nSym)-nFro(1:nSym)-nOccMax(:)-nDel(1:nSym)
  IndType(7,1:nSym) = nDel(1:nSym)
end if
if (nD == 1) then
  OrbName = 'SCFORB'
  if (KSDFT == 'SCF') then
    iWFtype = 2
  else
    Note = trim(Note)//' / '//trim(KSDFT)
    iWFtype = 3
  end if
  call WrVec_(OrbName,LuOut,What,nD-1,nSym,nBas,nBas,CMO(:,1),Dummy,OccNo(:,1),Dummy,E_Or(:,1),Dummy,IndType,Note,iWFtype)
# ifdef _HDF5_
  nZero = 0
  call mma_allocate(typestring,nnB)
  nSSh(1:nSym) = nBas(1:nSym)-nFro(1:nSym)-nOcc(1:nSym,1)-nDel(1:nSym)
  call orb2tpstr(NSYM,NBAS,NFRO,NOCC(1,1),NZERO,NZERO,NZERO,NSSH,NDEL,typestring)
  call mh5_put_dset(wfn_tpidx,typestring)
  call mma_deallocate(typestring)
  call mh5_put_dset(wfn_mocoef,CMO(1,1))
  call mh5_put_dset(wfn_occnum,OccNo(1,1))
  call mh5_put_dset(wfn_orbene,EOrb(1,1))
# endif
else
  OrbName = 'UHFORB'
  if (KSDFT == 'SCF') then
    iWFtype = 4
  else
    Note = trim(Note)//' / '//trim(KSDFT)
    iWFtype = 5
  end if
  call WrVec_(OrbName,LuOut,What,nD-1,nSym,nBas,nBas,CMO(:,1),CMO(:,2),OccNo(:,1),OccNo(:,2),E_Or(:,1),E_Or(:,2),IndType,Note, &
              iWFtype)
# ifdef _HDF5_
  nZero = 0
  call mma_allocate(typestring,nnB)
  nSSh(1:nSym) = nBas(1:nSym)-nFro(1:nSym)-nOcc(1:nSym,1)-nDel(1:nSym)
  call orb2tpstr(NSYM,NBAS,NFRO,NOCC(1,1),NZERO,NZERO,NZERO,NSSH,NDEL,typestring)
  call mh5_put_dset(wfn_tpidx_a,typestring)
  nSSh(1:nSym) = nBas(1:nSym)-nFro(1:nSym)-nOcc(1:nSym,2)-nDel(1:nSym)
  call orb2tpstr(NSYM,NBAS,NFRO,NOCC(1,2),NZERO,NZERO,NZERO,NSSH,NDEL,typestring)
  call mh5_put_dset(wfn_tpidx_b,typestring)
  call mh5_put_dset(wfn_mocoef_a,CMO(1,1))
  call mh5_put_dset(wfn_occnum_a,OccNo(1,1))
  call mh5_put_dset(wfn_orbene_a,EOrb(1,1))
  call mh5_put_dset(wfn_mocoef_b,CMO(1,2))
  call mh5_put_dset(wfn_occnum_b,OccNo(1,2))
  call mh5_put_dset(wfn_orbene_b,EOrb(1,2))
# endif
  iBas = 0
  IndType(1,1:nSym) = nFro(1:nSym)
  IndType(2,1:nSym) = 0
  IndType(3,1:nSym) = 0
  IndType(4,1:nSym) = 0
  IndType(5,1:nSym) = 0
  IndType(6,1:nSym) = nOrb(1:nSym)-nFro(1:nSym)-nDel(1:nSym)
  IndType(7,1:nSym) = nDel(1:nSym)
  do iSym=1,nSym
    do kBas=1,nBas(iSym)
      iBas = iBas+1
      if (Etan(iBas) > 1.99_wp) then
        IndType(2,iSym) = IndType(2,iSym)+1
        IndType(6,iSym) = IndType(6,iSym)-1
      else if (Etan(iBas) > 0.01_wp) then
        IndType(4,iSym) = IndType(4,iSym)+1
        IndType(6,iSym) = IndType(6,iSym)-1
      end if
    end do
  end do
  OrbName = 'UNAORB'
  Note = 'UHF natural orbitals'
  if (KSDFT == 'SCF') then
    iWFtype = 6
  else
    iWFtype = 7
  end if
  call WrVec_(OrbName,LuOut,What,0,nSym,nBas,nBas,CMOn,Dummy,Etan,Dummy,Epsn,Dummy,IndType,Note,iWFtype)
# ifdef _HDF5_
  IndTypeT = transpose(IndType)
  call orb2tpstr(NSYM,NBAS,NFRO,IndTypeT(:,2),NZERO,IndTypeT(:,4),NZERO,IndTypeT(:,6),NDEL,typestring)
  call mh5_put_dset(wfn_tpidx,typestring)
  call mma_deallocate(typestring)
  call mh5_put_dset(wfn_mocoef,CMOn)
  call mh5_put_dset(wfn_occnum,Etan)
  call mh5_put_dset(wfn_orbene,Epsn)
# endif
  call mma_deallocate(Epsn)
  call mma_deallocate(Etan)
  call mma_deallocate(CMOn)
end if

call mma_deallocate(E_Or)

! release Buffers for semi-direct SCF
if (DSCF) call ClsBuf()
! release SEWARD
if (DSCF .or. RF_On() .or. Langevin_On() .or. PCM_On() .or. Do_OFemb .or. Do_Tw .or. Do_Addc .or. &
#   ifdef _EFP_
    EFP_On() .or. &
#   endif
    (KSDFT /= 'SCF')) call ClsSew()

call mma_deallocate(FMaux,safe='*')
call mma_deallocate(NDSD,safe='*')
call mma_deallocate(DSc,safe='*')
#ifdef _EFP_
if (EFP_On()) call EFP_ShutDown(EFP_Instance)
#endif

call CWTime(TCpu2,TWall2)
TimFld(15) = TimFld(15)+(TCpu2-TCpu1)
TotCpu = max(TCpu2,0.1_wp)

! Write out timing informations
if (.not. NoProp) then
  nFldP = nFld-1
else
  nFldP = nFld-2
end if
Frmt = '(2x,A)'
if (iStatPRN > 0) then
  write(u6,*)
  call CollapseOutput(1,'Statistics and timing')
  write(u6,'(3X,A)') '---------------------'
  write(u6,*)
  write(u6,Frmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,Frmt) '   Part of the program                              CPU    fraction'
  write(u6,Frmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  do iFld=1,nFldP
    if ((iFld == 3) .or. (iFld == 4)) cycle
    write(u6,'(2x,A45,2f10.2)') NamFld(iFld),TimFld(iFld),TimFld(iFld)/TotCpu
  end do
  write(u6,*)
  write(u6,'(2x,A45,2F10.2)') NamFld(nFld),TotCpu
  write(u6,Frmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  call CollapseOutput(0,'Statistics and timing')
  write(u6,*)
end if

end subroutine FinalSCF
