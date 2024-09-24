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

subroutine Final()
!***********************************************************************
!                                                                      *
!     purpose: perform final calculations                              *
!                                                                      *
!***********************************************************************

use SCF_Arrays, only: Dens, OneHam, Ovrlp, TwoHam, CMO, EOrb, FockAO, OccNo, KntE, MssVlc, Darwin
use InfSCF, only: nD
#ifdef _EFP_
use EFP_Module, only: EFP_Instance
use EFP, only: EFP_Shutdown
#endif
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use SCFWFn, only: wfn_mocoef, wfn_mocoef_a, wfn_mocoef_b, wfn_occnum, wfn_occnum_a, wfn_occnum_b, wfn_orbene, wfn_orbene_a, &
                  wfn_orbene_b, wfn_tpidx, wfn_tpidx_a, wfn_tpidx_b
#endif
use Interfaces_SCF, only: dOne_SCF
use OFembed, only: Do_OFemb, FMaux, NDSD
#ifdef _FDE_
use Embedding_Global, only: embPot, embWriteEsp
#endif
use SpinAV, only: DSc
use InfSCF, only: nBT, nDens, DMOMax, FMOMax, kIVO, MaxBas, nSym, KSDFT, EneV, Falcon, iPrint, NoProp, DSCF, TotCPU, nFld, &
                  iStatPrn, E1V, E2V, FThr, iPrForm, MaxBXO, Name, NamFld, nBas, nBB, nBO, nDel, nFro, nIter, nIterP, nnB, nnO, &
                  nOcc, nOrb, TimFld
#ifdef _FDE_
use InfSCF, only: nAtoms
#endif
use Constants, only: Zero, One, Two
use stdalloc, only: mma_allocate, mma_deallocate
use Files, only: LuOut
use AddCorr, only: Do_Addc, Do_Tw

implicit none
#ifdef _EFP_
logical, external :: EFP_On
#endif
! Define local variable
integer iD, iRC, iOpt, iSymLb, iFock, jFock, iCMO, iVirt, jVirt, ij, iBas, jBas, iSym, kl, lk, iRef, jRef, iiOrb, iOrb, nOccMax, &
        nOccMin, iWFType, kBas, iFld
real*8 TCPU1, TCPU2, Dummy, TWall1, TWall2
logical FstItr
character(len=8) RlxLbl, Method
character(len=60) Fmt
character(len=128) OrbName
logical RF_On, Langevin_On, PCM_On
character(len=80) Note
character(len=8) What
integer IndType(7,8)
real*8, dimension(:), allocatable :: Temp, CMOn, Etan, Epsn
real*8, dimension(:,:), allocatable :: GVFck, Scrt1, Scrt2, DMat, EOr
#ifdef _HDF5_
#include "Molcas.fh"
character(len=1), allocatable :: typestring(:)
integer nSSh(mxSym), nZero(mxSym)
integer i
integer IndTypeT(8,7)
#endif
integer nFldP
dimension Dummy(1)

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
  call dcopy_(nBT*nD,Dens(1,1,1),1,Dens(1,1,nDens),1)
  call dcopy_(nBT*nD,TwoHam(1,1,1),1,TwoHam(1,1,nDens),1)

  DMOMax = Zero
  FMOMax = Zero

else

  ! Compute improved virtuals (if needed)
  if (kIvo == 1) then
    do iD=1,nD
      call IvoGen(OneHam,nBT,CMO(1,iD),nBO,EOrb(1,iD),nnO,nOcc(1,iD))
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
  call FZero(Temp,nBT+4)

  do iD=1,nD
    call DCopy_(nBT,FockAO(1,iD),1,Temp,1)
    iRc = -1
    iOpt = 0
    RlxLbl = 'Fock Op '
    if (iD == 2) RlxLbl = 'Fock Op2'
    iSymLb = 1
    call WrOne(iRc,iOpt,RlxLbl,1,Temp,iSymLb)
    if (iRc /= 0) then
      write(6,*) 'Final: Error writing on ONEINT'
      write(6,'(A,A)') 'RlxLbl=',RlxLbl
      call Abend()
    end if
  end do
  call mma_deallocate(Temp)

  ! Construct the generalized variational Fock matrix

  call mma_allocate(GVFck,nBT,nD,Label='GVFck')
  call mma_allocate(Scrt1,MaxBas**2,nD,Label='Scrt1')
  call mma_allocate(Scrt2,MaxBxO,nD,Label='Scrt2')

  iFock = 1
  jFock = 1
  iCMo = 1
  do iSym=1,nSym
    if (nOrb(iSym) <= 0) cycle

    do iD=1,nD
      call Square(FockAO(jFock,iD),Scrt1(1,iD),1,nBas(iSym),nBas(iSym))
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
        call dcopy_(nOrb(iSym)-nOcc(iSym,iD),[Zero],0,Scrt1(jVirt,iD),1)
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

    iFock = iFock+nBas(iSym)*(nBas(iSym)+1)/2
    jFock = jFock+nBas(iSym)*(nBas(iSym)+1)/2
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
      call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(1,1),nBB,OccNo(1,1),DMat(1,1),.true.)
    else
      call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(1,1),nBB,OccNo(1,1),DMat(1,1),.true.)
      call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(1,2),nBB,OccNo(1,2),DMat(1,2),.false.)
      call Put_dArray('D1aoVar',Dummy,0) ! Undefined the field.
    end if

    call DensAB(nBT,1,nD,DMat)
    call mma_deallocate(DMat)
  end if

  if (nD == 2) then
    call mma_allocate(CMOn,nBB,Label='CMOn')
    call mma_allocate(Etan,nnB,Label='Etan')
    call mma_allocate(Epsn,nnB,Label='Epsn')
    call NatoUHF(Dens(1,1,1),Dens(1,2,1),FockAO(1,1),FockAO(1,2),nBT,CMO(1,1),nBB,Ovrlp,CMOn,Etan,Epsn,nnB,nSym,nBas,nOrb)
    call PadCMO(CMOn,CMOn,nSym,nBas,nOrb)
    call PadEor(Etan,Etan,nSym,nBas,nOrb)
    call PadEor(Epsn,Epsn,nSym,nBas,nOrb)
  end if

  ! Add generalized Fock matrix
  if (nD == 2) then
    call DaXpY_(nBT,One,GVFck(1,2),1,GVFck(1,1),1)
  else
    call DScal_(nBT,Two,GvFck(1,1),1)
  end if
  call Put_dArray('FockOcc',GVFck(1,1),nBT)
  call mma_deallocate(GVFck)

  ! Add SCF orbital energies
  call Put_dArray('OrbE',EOrb(1,1),nnB)
  if (nD == 2) call Put_dArray('OrbE_ab',EOrb(1,2),nnB)

  call Put_dScalar('Thrs',Fthr)
end if

#ifdef _FDE_
! Embedding
if (embPot .and. embWriteEsp) call embPotOutput(nAtoms,Dens)
#endif

! t.t.;
! store Fock matrix in Runfile in case of fragment calculation
if (Falcon) then
  !write(6,*) 'Fock matrix is written in RunFile.'
  !write(6,*) 'fck:'
  !write(6,*) 'nbt=',nbt
  !write(6,*) 'ndens=',ndens
  !write(6,*) (FockAO(itt),itt=1,nbt)
  call Put_dArray('Fragment_Fock',FockAO(1,1),nBT)
end if

! t.t.; end
! Final printout
! original PrFin was splitted to 3 parts to run in UHF mode

do iSym=1,nSym
  nDel(iSym) = nBas(iSym)-nOrb(iSym)
end do
call PrFin0(Dens(1,1,nDens),Dens(1,2,nDens),nBT,EOrb(1,1),nnB,CMO(1,1),nBO,KntE)

do iD=1,nD
  call PrFin(OneHam,Ovrlp,Dens(1,iD,nDens),TwoHam(1,iD,nDens),nBT,EOrb(1,iD),OccNo(1,iD),nnB,CMO(1,iD),nBO,Note,iD-1,MssVlc,Darwin)

  call PrFin2(Ovrlp,nBT,OccNo(1,iD),nnB,CMO(1,iD),nBO,Note)
end do

if ((iPrint >= 2) .and. (nD == 2)) &
  call PriMO('Natural orbitals',.true.,.true.,Zero,2.0d9,nSym,nBas,nOrb,Name,Epsn,Etan,CMOn,iPrForm)
! Calculate expectation values
if (.not. NoProp) then
  if (iPrint >= 3) write(6,'(/6X,A)') 'Expectation values of various operators'
  call Prpt()
end if

! make a fix for energies for deleted orbitals
call mma_allocate(EOr,nnB,nD,Label='EOr')

iRef = 1
jRef = 1
do iSym=1,nSym
  iiOrb = nOrb(iSym)-nDel(iSym)
  do iOrb=1,iiOrb
    do iD=1,nD
      EOr(iRef,iD) = EOrb(jRef,iD)
    end do
    iRef = iRef+1
    jRef = jRef+1
  end do
  do iOrb=1,nDel(iSym)
    do iD=1,nD
      EOr(iRef,iD) = 1000
    end do
    iRef = iRef+1
  end do
end do

if (nD == 1) then
  do iSym=1,nSym
    IndType(1,iSym) = nFro(iSym)
    IndType(2,iSym) = nOcc(iSym,1)
    IndType(3,iSym) = 0
    IndType(4,iSym) = 0
    IndType(5,iSym) = 0
    IndType(6,iSym) = nOrb(iSym)-nFro(iSym)-nOcc(iSym,1)-nDel(iSym)
    IndType(7,iSym) = nDel(iSym)
  end do
else
  do iSym=1,nSym
    nOccMax = max(nOcc(iSym,1),nOcc(iSym,2))
    nOccMin = min(nOcc(iSym,1),nOcc(iSym,2))
    IndType(1,iSym) = nFro(iSym)
    IndType(2,iSym) = nOccMin
    IndType(3,iSym) = 0
    IndType(4,iSym) = nOccMax-nOccMin
    IndType(5,iSym) = 0
    IndType(6,iSym) = nOrb(iSym)-nFro(iSym)-nOccMax-nDel(iSym)
    IndType(7,iSym) = nDel(iSym)
  end do
end if
if (nD == 1) then
  OrbName = 'SCFORB'
  if (KSDFT == 'SCF') then
    iWFtype = 2
  else
    Note = trim(Note)//' / '//trim(KSDFT)
    iWFtype = 3
  end if
  call WrVec_(OrbName,LuOut,What,nD-1,nSym,nBas,nBas,CMO(1,1),Dummy,OccNo(1,1),Dummy,EOr(1,1),Dummy,IndType,Note,iWFtype)
# ifdef _HDF5_
  nZero = 0
  call mma_allocate(typestring,nnB)
  do i=1,nSym
    nSSh(i) = nBas(i)-nFro(i)-nOcc(i,1)-nDel(i)
  end do
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
  call WrVec_(OrbName,LuOut,What,nD-1,nSym,nBas,nBas,CMO(1,1),CMO(1,2),OccNo(1,1),OccNo(1,2),EOr(1,1),EOr(1,2),IndType,Note,iWFtype)
# ifdef _HDF5_
  nZero = 0
  call mma_allocate(typestring,nnB)
  do i=1,nSym
    nSSh(i) = nBas(i)-nFro(i)-nOcc(i,1)-nDel(i)
  end do
  call orb2tpstr(NSYM,NBAS,NFRO,NOCC(1,1),NZERO,NZERO,NZERO,NSSH,NDEL,typestring)
  call mh5_put_dset(wfn_tpidx_a,typestring)
  do i=1,nSym
    nSSh(i) = nBas(i)-nFro(i)-nOcc(i,2)-nDel(i)
  end do
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
  do iSym=1,nSym
    IndType(1,iSym) = nFro(iSym)
    IndType(2,iSym) = 0
    IndType(3,iSym) = 0
    IndType(4,iSym) = 0
    IndType(5,iSym) = 0
    IndType(6,iSym) = nOrb(iSym)-nFro(iSym)-nDel(iSym)
    IndType(7,iSym) = nDel(iSym)
    do kBas=1,nBas(iSym)
      iBas = iBas+1
      if (Etan(iBas) > 1.99d0) then
        IndType(2,iSym) = IndType(2,iSym)+1
        IndType(6,iSym) = IndType(6,iSym)-1
      else if (Etan(iBas) > 0.01d0) then
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

call mma_deallocate(EOr)

! release Buffers for semi-direct SCF
if (DSCF) call ClsBuf()
! release SEWARD
if (DSCF .or. RF_On() .or. Langevin_On() .or. PCM_On() .or. Do_OFemb .or. Do_Tw .or. Do_Addc .or. &
#   ifdef _EFP_
    EFP_On() .or. &
#   endif
    (KSDFT /= 'SCF')) call ClsSew()

if (allocated(FMaux)) call mma_deallocate(FMaux)
if (allocated(NDSD)) call mma_deallocate(NDSD)
if (allocated(DSc)) call mma_deallocate(DSc)
#ifdef _EFP_
if (EFP_On()) call EFP_ShutDown(EFP_Instance)
#endif

call CWTime(TCpu2,TWall2)
TimFld(15) = TimFld(15)+(TCpu2-TCpu1)
TotCpu = max(TCpu2,0.1d0)

! Write out timing informations
if (.not. NoProp) then
  nFldP = nFld-1
else
  nFldP = nFld-2
end if
Fmt = '(2x,A)'
if (iStatPRN > 0) then
  write(6,*)
  call CollapseOutput(1,'Statistics and timing')
  write(6,'(3X,A)') '---------------------'
  write(6,*)
  write(6,Fmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  write(6,Fmt) '   Part of the program                              CPU    fraction'
  write(6,Fmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  do iFld=1,nFldP
    if ((iFld == 11) .or. (iFld == 12)) cycle
    write(6,'(2x,A45,2f10.2)') NamFld(iFld),TimFld(iFld),TimFld(iFld)/TotCpu
  end do
  write(6,*)
  write(6,'(2x,A45,2F10.2)') NamFld(nFld),TotCpu
  write(6,Fmt) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  call CollapseOutput(0,'Statistics and timing')
  write(6,*)
end if

end subroutine Final
