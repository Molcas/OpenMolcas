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
! Copyright (C) Ben Swerts                                             *
!               2016, Liviu Ungur                                      *
!***********************************************************************

subroutine Drv2El_FAIEMP()
!***********************************************************************
!                                                                      *
!  Object: driver for the central-fragment interaction 2-electron      *
!          integrals (based on drv2el_3center_RI and drv2el_scf)       *
!                                                                      *
!     Author: Ben Swerts                                               *
!   Modified: Liviu Ungur                                              *
!***********************************************************************

#ifdef _DEBUGPRINT_
use k2_arrays, only: pDq, pFq
#endif
use Basis_Info, only: dbsc, nBas, nBas_Frag, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iOper
use Gateway_Info, only: ThrInt, CutInt
use Int_Options, only: Disc, Disc_Mx, DoFock, DoIntegrals, ExFac, FckNoClmb, FckNoExch, PreSch, Thize, W2Disc
use Integral_Interfaces, only: DeDe_SCF, Int_PostProcess, No_Routine
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: nTInt = 1
integer(kind=iwp) :: i, iCnt, iCnttp, iComp, iDpos, iFpos, iIrrep, ijS, iOpt, iRC, iS, j, jS, klS, kS, lOper, lS, maxDens, mdc, &
                     mDens, nBas_Valence(0:7), nBasC, nBT, nBVT, nBVTi, nFock, nij, nOneHam, nSkal, nSkal_Valence
real(kind=wp) :: A_int, Dtst, P_Eff, TCpu1, TCpu2, ThrAO, TInt(nTInt), TMax_all, TWall1, TWall2
logical(kind=iwp) :: DoGrad, EnergyWeight, Indexation, lNoSkip
character(len=8) :: Label
integer(kind=iwp), allocatable :: ij(:)
real(kind=wp), allocatable, target :: Dens(:), Fock(:)
real(kind=wp), allocatable :: DMax(:,:), FragDensSO(:), OneHam(:), TMax(:,:)
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iFD, nFD
character(len=80) :: Line
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
call xFlush(u6)
! set variables in module Int_Options
DoIntegrals = .false.
DoFock = .true.
FckNoClmb = .false. ! Default value
FckNoExch = .false. ! Default value
ExFac = One
Thize = 1.0e-6_wp
W2Disc = .true.
PreSch = .false.
Disc_Mx = Zero      ! Default value
Disc = Zero         ! Default value
Int_PostProcess => No_Routine

! Handle both the valence and the fragment basis set

call Set_Basis_Mode('Valence')
call Nr_Shells(nSkal_Valence)
call Free_iSD()
call Set_Basis_Mode('WithFragments')
call SetUp_iSD()
nBT = 0
nBVT = 0
do i=0,nIrrep-1
  nBas_Valence(i) = nBas(i)
  nBVT = nBVT+nBas(i)*(nBas(i)+1)/2
  nBas(i) = nBas(i)+nBas_Frag(i)
  nBT = nBT+nBas(i)*(nBas(i)+1)/2
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!---  Construct custom density matrix

call mma_allocate(Dens,nBT,Label='Dens')
call mma_allocate(Fock,nBT,Label='Fock')
! Valence part is zero
Dens(:) = Zero
Fock(:) = Zero
! Each fragment needs it's (symmetrized) density matrix added along the
! diagonal.
! This density matrix first has to be constructed from the MO coeffs
! so allocate space for the largest possible density matrix
maxDens = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%nFragType > 0) maxDens = max(maxDens,dbsc(iCnttp)%nFragDens*(dbsc(iCnttp)%nFragDens+1)/2)
end do
call mma_allocate(FragDensSO,maxDens,label='FragDSO')

iDpos = 1 ! position in the total density matrix
do iIrrep=0,nIrrep-1
  nBasC = nBas_Valence(iIrrep)
  iDpos = iDpos+nBasC*(nBasC+1)/2
  mdc = 0
  do iCnttp=1,nCnttp
    if (dbsc(iCnttp)%nFragType <= 0) then
      mdc = mdc+dbsc(iCnttp)%nCntr
    else

      ! construct the density matrix
      EnergyWeight = .false.
      !call MakeDens(dbsc(iCnttp)%nFragDens,dbsc(iCnttp)%nFragEner,dbsc(iCnttp)%FragCoef,dbsc(iCnttp)%FragEner,EnergyWeight, &
      !              FragDensAO)
      call MakeDens(dbsc(iCnttp)%nFragDens,dbsc(iCnttp)%nFragEner,dbsc(iCnttp)%FragCoef,dbsc(iCnttp)%FragEner,EnergyWeight, &
                    FragDensSO)
      ! create the symmetry adapted version if necessary
      ! (fragment densities are always calculated without symmetry)
#     ifdef _DEBUGPRINT_
      call TriPrt('Fragment density',' ',FragDensSO,dbsc(iCnttp)%nFragDens)
#     endif

      do iCnt=1,dbsc(iCnttp)%nCntr
        mdc = mdc+1
        ! only add fragment densities that are active in this irrep
        ! => the following procedure still has to be verified thoroughly
        !    but appears to be working
        if (iand(dc(mdc)%iChCnt,iIrrep) == iOper(iIrrep)) then
          ! add it at the correct location in the large custom density matrix
          iFpos = 1
          ! position in fragment density matrix
          do i=1,dbsc(iCnttp)%nFragDens
            iDpos = iDpos+nBasC
            do j=0,i-1
              Dens(iDpos+j) = FragDensSO(iFpos+j)
            end do
            iDpos = iDpos+i
            iFpos = iFpos+i
          end do
          nBasC = nBasC+dbsc(iCnttp)%nFragDens
        end if
      end do
    end if
  end do
end do
#ifdef _DEBUGPRINT_
iFD = 0
do iIrrep=0,nIrrep-1
  nFD = nBas(iIrrep)*(nBas(iIrrep)+1)/2
  call TriPrt('Combined density',' ',Dens(iFD+1:iFD+nFD),nBas(iIrrep))
  iFD = iFD+nFD
end do
#endif
call mma_deallocate(FragDensSO)
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Desymmetrize the custom density matrix
!
! Should be possible to reduce the storage space for the Fock matrix
! and only store the top nBas_Valence(0) rows (non-symmetry case tested)

call AlloK2()
call DeDe_SCF(Dens,Fock,nBT,mDens)
#ifdef _DEBUGPRINT_
if (nIrrep == 1) then
  call RecPrt('Desymmetrized Density:',' ',pDq,nBas(0),nBas(0))
else
  iFD = 0
  do iIrrep=0,nIrrep-1
    nFD = nBas(iIrrep)*(nBas(iIrrep)+1)/2
    call TriPrt('Desymmetrized density',' ',pDq(iFD+1:iFD+nFD),nBas(iIrrep))
    iFD = iFD+nFD
  end do
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation. Do not generate
! tables for indexation.

ThrAO = Zero
Indexation = .false.
DoGrad = .false.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!                                                                      *
!***********************************************************************
!                                                                      *

ThrInt = CutInt   ! Integral neglect threshold from SCF
!                                                                      *
!***********************************************************************
!                                                                      *
!---  Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal,label='TMax')
call Shell_MxSchwz(nSkal,TMax)
TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do
call mma_allocate(DMax,nSkal,nSkal,label='DMax')
call Shell_MxDens(Dens,DMax,nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(ij,nSkal*(nSkal+1),label='ij')
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      ij(nij*2-1) = iS
      ij(nij*2) = jS
    end if
  end do
end do
P_Eff = real(nij,kind=wp)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)

! Now do a quadruple loop over shells

!ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
ijS = 1
iS = ij(ijS*2-1)
jS = ij(ijS*2)
!klS = int(TskLw-real(ijS,kind=wp)*(real(ijS,kind=wp)-One)/Two)
klS = 1
kS = ij(klS*2-1)
lS = ij(klS*2)
do
  if (ijS > int(P_Eff)) exit
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! density prescreening (results in iS > nSkal_Valence)
  A_int = TMax(iS,jS)*TMax(kS,lS)
  Dtst = max(DMax(is,ls)*Quart,DMax(is,ks)*Quart,DMax(js,ls)*Quart,DMax(js,ks)*Quart,DMax(is,js),DMax(ks,ls))
  lNoSkip = A_int*Dtst >= ThrInt
  ! only calculate needed integrals and only update the valence part of the
  ! Fock matrix (iS > nSkal_Valence, lS <= nSkal_Valence, jS and kS
  ! belonging to different regions)
  if (jS <= nSkal_Valence) then
    lNoSkip = lNoSkip .and. (kS > nSkal_Valence)
  else
    lNoSkip = lNoSkip .and. (kS <= nSkal_Valence)
  end if
  lNoSkip = lNoSkip .and. (lS <= nSkal_Valence)

  if (lNoSkip) then
    call Eval_IJKL(iS,jS,kS,lS,TInt,nTInt)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Drv2El_FAIEMP: for iS, jS, kS, lS =',is,js,ks,ls
    if (nIrrep == 1) then
      call RecPrt('updated Fock',' ',pFq,nBas(0),nBas(0))
    else
      iFD = 0
      do iIrrep=0,nIrrep-1
        nFD = nBas(iIrrep)*(nBas(iIrrep)+1)/2
        call TriPrt('updated Fock',' ',pFq(iFD+1:iFD+nFD),nBas(iIrrep))
        iFD = iFD+nFD
      end do
    end if
#   endif
  end if
  klS = klS+1
  if (klS > ijS) then
    ijS = ijS+1
    klS = 1
  end if
  iS = ij(ijS*2-1)
  jS = ij(ijS*2)
  kS = ij(klS*2-1)
  lS = ij(klS*2)
end do

call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(ij)
call mma_deallocate(DMax)
call mma_deallocate(TMax)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Terminate integral environment.

call Term_Ints()

call Free_DeDe(Dens,Fock,nBT)

call mma_deallocate(Dens)
Int_PostProcess => null()
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*)
write(u6,'(a)') 'SO Integrals of type Frag2El Component 1'
iFD = 0
do iIrrep=0,nIrrep-1
  write(Line,'(1X,A,I1)') ' Diagonal Symmetry Block ',iIrrep+1
  nFD = nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
  call TriPrt(Line,' ',Fock(iFD+1:iFD+nFD),nBas_Valence(iIrrep))
  iFD = iFD+nFD
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Write the results to the one electron integral file
!
! read the one electron hamiltonian
Label = 'OneHam  '
iRC = -1
iOpt = 0
iComp = 1
call mma_allocate(OneHam,nBVT+4,label='OneHam')
call RdOne(iRC,iOpt,Label,iComp,OneHam,lOper)
if (iRC /= 0) then
  write(u6,*) 'Drv2El_FAIEMP: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
! add the calculated results
nOneHam = 1 ! counter in the OneHam matrices (small)
nFock = 1   ! counter in the Fock matrices (larger)
do iIrrep=0,nIrrep-1
  nBVTi = nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
  call daxpy_(nBVTi,One,Fock(nFock),1,OneHam(nOneHam),1)
  nOneHam = nOneHam+nBVTi
  nFock = nFock+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

! write out the results
#ifdef _DEBUGPRINT_
iFD = 0
do iIrrep=0,nIrrep-1
  nFD = nBas_Valence(iIrrep)*(nBas_Valence(iIrrep)+1)/2
  call TriPrt('OneHam at end',' ',OneHam(iFD+1:iFD+nFD),nBas_Valence(iIrrep))
  iFD = iFD+nFD
end do
#endif
iRC = -1
call WrOne(iRC,iOpt,Label,1,OneHam,lOper)
if (iRC /= 0) then
  write(u6,*) 'Drv2El_FAIEMP: Error writing to ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
iRC = -1
Label = 'OneHam 0'
call WrOne(iRC,iOpt,Label,1,OneHam,lOper)
if (iRC /= 0) then
  write(u6,*) 'Drv2El_FAIEMP: Error writing to ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

! cleanup
call mma_deallocate(OneHam)
call mma_deallocate(Fock)
!                                                                      *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
call Free_iSD()
call Set_Basis_Mode('Valence')
call SetUp_iSD()
do i=0,nIrrep-1
  nBas(i) = nBas_Valence(i)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call Init_Int_Options()
return

end subroutine Drv2El_FAIEMP
