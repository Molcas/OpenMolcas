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

subroutine Gamma_new(Int1,Int2,Int1_2,Int2_2,Scr1)

#include "intent.fh"

use MBPT2_Global, only: CMO, EOcc, EVir, mAdOcc, mAdVir, nBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
#include "corbinf.fh"
real(kind=wp), intent(_OUT_) :: Int1(*), Int2(*), Int1_2(*), Int2_2(*), Scr1(*)
integer(kind=iwp) :: i, iA, iAdrBin, iAdrGam, iAdrRdBin, iB, iBB, iBin, iBinLength, iBinSize, iBlock, iI, iIA, iJ, iJJ, iKap, &
                     iLam, iLamKap1, iLamKap2, iLastAdr, iLen, iMaxBas, iMaxBasProd, iMaxOccVir, iMemAvail, iMemNeeded, iMu, &
                     iMuNu1, iMuNu2, iNextAdr, iNextX, iNu, iOff, iOffCMO(nSym), iOffCMO_o(nSym), iOffCMO_v(nSym), iRec, iSize, &
                     iSym, iSym1, iSym2, iSym_A, iSym_B, iSym_C, iSym_D, iTriMuNu, iType, lAllBins, lCMO_o, lCMO_v, lMax, LuBin, &
                     LuGam, nA, nA2, nABCD, nB, nB2, nBins, nBlocks, nI, nI2, nJ, nJ2, nLam, nNO, nNu, nNV, nTOrb(nSym)
real(kind=wp) :: EDenom, Tiajb, xiajb, xibja
logical(kind=iwp) :: LoadZeros, NonZeroSym(4), Triangular
integer(kind=iwp), allocatable :: iTable(:,:)
real(kind=wp), allocatable :: Bin(:,:,:), Bin2(:,:,:), CMO_o(:), CMO_v(:), Temp1(:), Temp1_2(:), Temp2(:), Temp2_2(:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
iAdrGam = 1
LoadZeros = .false.

iMaxBas = 0
do iSym=1,nSym
  nTOrb(iSym) = nOrb(iSym)+nDel(iSym)
  iMaxBas = max(iMaxBas,nBas(iSym))
end do

! Setup one CMO for occupied and one for virtual orbitals together
! with a corresponding symmetryoffset.

iOffCMO(1) = 0
iOffCMO_o(1) = 0
iOffCMO_v(1) = 0
do iSym=2,nSym
  iOffCMO(iSym) = iOffCMO(iSym-1)+nTOrb(iSym-1)*nBas(iSym-1)
  iOffCMO_o(iSym) = iOffCMO_o(iSym-1)+nOcc(iSym-1)*nBas(iSym-1)
  iOffCMO_v(iSym) = iOffCMO_v(iSym-1)+nExt(iSym-1)*nBas(iSym-1)
end do

! Calculate total length of CMO-matrices in all irreps.

lCMO_o = 0
lCMO_v = 0
do iSym=1,nSym
  lCMO_o = lCMO_o+nOcc(iSym)*nBas(iSym)
  lCMO_v = lCMO_v+nExt(iSym)*nBas(iSym)
end do

! Allocate memory for a virtual and an occupied CMO-matrix

call mma_allocate(CMO_o,lCMO_o,label='CMO_o')
call mma_allocate(CMO_v,lCMO_v,label='CMO_v')

! Copy CMO to CMO_o and CMO_v

do iSym=1,nSym
  iOff = nBas(iSym)*nFro(iSym)+1
  nNO = nBas(iSym)*nOcc(iSym)
  nNV = nBas(iSym)*nExt(iSym)

  call dCopy_(nNO,CMO(iOffCMO(iSym)+iOff),1,CMO_o(iOffCMO_o(iSym)+1),1)

  iOff = iOff+nNO

  call dCopy_(nNV,CMO(iOffCMO(iSym)+iOff),1,CMO_v(iOffCMO_v(iSym)+1),1)
end do

#ifdef _DEBUGPRINT_
! Print the elements of the Full CMO-matrices as well as CMO_o and CMO_v.
do iSym=1,nSym
  call RecPrt('Full CMO',' ',CMO(iOffCMO(iSym)+1),nBas(iSym),nTOrb(iSym))
end do
do iSym=1,nSym
  call RecPrt('Occupied CMO',' ',CMO_o(iOffCMO_o(iSym)+1),nBas(iSym),nOcc(iSym))
end do
do iSym=1,nSym
  call RecPrt('Virtual CMO',' ',CMO_v(iOffCMO_v(iSym)+1),nBas(iSym),nExt(iSym))
end do
#endif

! Setup indices for symmetryblocks to make it compatible with
! the aces-routine used for reading the gammas later.

if (nSym == 8) nBlocks = 106
if (nSym == 4) nBlocks = 19
if (nSym == 2) nBlocks = 4
if (nSym == 1) nBlocks = 1

! Construct a table for block information.

call mma_allocate(iTable,6,nBlocks,label='iTable')
call Gamma_Blocks(iTable,nBlocks,nSym)

! Open a file for writing gammas.

LuGam = IsFreeUnit(10)
call DaName_MF_WA(LuGam,'LuGam')

! Open a file to store bins of half-transformed gammas.

LuBin = IsFreeUnit(11)
call DaName_MF_WA(LuBin,'TmpBin')

! SETUP OF AMPLITUDE BATCHES
! (The handling is somewhat simple when nBlocks=1 since we only have
!  blocks with one symmetry)
!
! Check how large the occ*vir-products and bas**2-products may be.

iMaxBasProd = 0
iMaxOccVir = 0
do iSym1=1,nSym
  do iSym2=1,iSym1
    iMaxBasProd = max(iMaxBasProd,nBas(iSym1)*nBas(iSym2))
    iMaxOccVir = max(iMaxOccVir,nOcc(iSym1)*nExt(iSym2))
    iMaxOccVir = max(iMaxOccVir,nExt(iSym1)*nOcc(iSym2))
  end do
end do

! Setup two pairs of temporary vectors to use for transformation
! and symmetrization

call mma_allocate(Temp1,iMaxBasProd,label='Temp1')
call mma_allocate(Temp2,iMaxBasProd,label='Temp2')
if (nBlocks /= 1) then
  call mma_allocate(Temp1_2,iMaxBasProd,label='Temp1_2')
  call mma_allocate(Temp2_2,iMaxBasProd,label='Temp2_2')
end if

! Check max available memory

call mma_maxDBLE(lMax)

! Calculate how much space is needed to store all full bins
! on disk.

if (nBlocks == 1) then
  iMemNeeded = 2*iMaxBasProd*(iMaxOccVir+1)
else
  iMemNeeded = 4*iMaxBasProd*(iMaxOccVir+1)
end if

! Take maximum one third of the memory for each bin-complex
! and leave one third for other tasks.

if (nBlocks == 1) then
  iMemAvail = 2*lMax/3
else
  iMemAvail = lMax/3
end if

! Make the bins as large as possible but not larger than is
! needed to keep all bins in memory.

lAllBins = min(iMemAvail,iMemNeeded)
iBinLength = lAllBins/iMaxBasProd

! Allocate two bins

!FIXME: Mixed types in Bin/Bin2
call mma_allocate(Bin,[1,2],[0,iBinLength/2-1],[1,iMaxBasProd],label='Bins1')
if (nBlocks /= 1) then
  call mma_allocate(Bin2,[1,2],[0,iBinLength/2-1],[1,iMaxBasProd],label='Bins2')
end if

! The construction and transformation of Tiajb are now made
! looping over symmetryblocks

do iBlock=1,nBlocks
  iType = iTable(1,iBlock)

  if ((iType == 1) .or. (iType == 2)) then
    Triangular = .true.
  else
    Triangular = .false.
  end if

  ! nBas, nOcc etc. are defined as nBas(1:8) here in MBPT2
  ! and as nBas(0:7) in integral_util and alaska. (which might not
  ! be optimal...)
  ! The weird order here is because I accidently constructed
  ! (kap lam|mu nu) when I needed (mu nu|kap lam), with this order
  ! it will be compatible with src/integral_util/read_blocks.f
  iSym_A = iTable(4,iBlock)+1
  iSym_B = iTable(5,iBlock)+1
  iSym_C = iTable(2,iBlock)+1
  iSym_D = iTable(3,iBlock)+1

  ! Number of orbitals in (ia|jb)
  nI = nOcc(iSym_A)
  nA = nExt(iSym_B)
  nJ = nOcc(iSym_C)
  nB = nExt(iSym_D)

  ! Number of orbitals in (ai|bj) (nX and nX1 can be
  ! combined to get (ia|bj) and (ai|jb) which are needed for
  ! full symmetrization)
  nI2 = nOcc(iSym_B)
  nA2 = nExt(iSym_A)
  nJ2 = nOcc(iSym_D)
  nB2 = nExt(iSym_C)

  ! Initialize adress for storing Bins on disk.
  iAdrBin = 1
  ! Setup the number of bins needed.
  if (Triangular) then
    nBins = nBas(iSym_C)*(nBas(iSym_C)+1)/2
  else
    nBins = nBas(iSym_C)*nBas(iSym_D)
  end if

  ! Initialize the bins to have length 0 and adress -1 to
  ! next element
  Bin(1,0,1:nBins) = Zero
  Bin(2,0,1:nBins) = -One
  if (nBlocks /= 1) then
    Bin2(1,0,1:nBins) = Zero
    Bin2(2,0,1:nBins) = -One
  end if

  ! Check which symmetry combinations that are present.
  NonZeroSym(:) = .true.
  if (nI*nA*nJ*nB == 0) NonZeroSym(1) = .false.
  if (nI*nA*nJ2*nB2 == 0) NonZeroSym(2) = .false.
  if (nI2*nA2*nJ*nB == 0) NonZeroSym(3) = .false.
  if (nI2*nA2*nJ2*nB2 == 0) NonZeroSym(4) = .false.

  ! Skip cases where nABCD is zero. If nABCD is nonzero but
  ! nIJAB (different A and B) is zero we jump to the end and load
  ! some zeros.

  if (Triangular) then
    nABCD = (nBas(iSym_A)*(nBas(iSym_A)+1)/2)*(nBas(iSym_C)*(nBas(iSym_C)+1)/2)
  else
    nABCD = nBas(iSym_A)*nBas(iSym_B)*nBas(iSym_C)*nBas(iSym_D)
  end if

  if (nABCD == 0) cycle
  if (((.not. NonZeroSym(1)) .and. Triangular) .or. &
      (.not. (NonZeroSym(1) .or. NonZeroSym(2) .or. NonZeroSym(3) .or. NonZeroSym(4)))) then
    LoadZeros = .true.
  else

    ! Start the loop to construct Tiajb.
    do iI=1,nI
      do iA=1,nA
        if (NonZeroSym(1)) then
          call Exch(iSym_D,iSym_A,iSym_C,iSym_B,iI+nFro(iSym_A),iA+nFro(iSym_B)+nOcc(iSym_B),Int1,Scr1)
          call Coul(iSym_D,iSym_C,iSym_A,iSym_B,iI+nFro(iSym_A),iA+nFro(iSym_B)+nOcc(iSym_B),Int2,Scr1)
#         ifdef _DEBUGPRINT_
          write(u6,*) ' *  I,A = ',iI,iA
          call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSym_D)+nDel(iSym_D),nOrb(iSym_C)+nDel(iSym_C))
          write(u6,*) ' *  I,A = ',iI,iA
          call RecPrt('Int2:','(8F10.6)',Int2,nOrb(iSym_D)+nDel(iSym_D),nOrb(iSym_C)+nDel(iSym_C))
#         endif
        end if

        ! When Symmetry of I and A differs we need to load both
        ! possible (IA|xx) to have (mu nu| xx) and (nu mu|xx)
        ! available at the same time for symmetrization.
        if (NonZeroSym(2) .and. (.not. Triangular)) then
          call Exch(iSym_C,iSym_A,iSym_D,iSym_B,iI+nFro(iSym_A),iA+nFro(iSym_B)+nOcc(iSym_B),Int1_2,Scr1)
          call Coul(iSym_C,iSym_D,iSym_A,iSym_B,iI+nFro(iSym_A),iA+nFro(iSym_B)+nOcc(iSym_B),Int2_2,Scr1)
#         ifdef _DEBUGPRINT_
          write(u6,*) ' *  I,A = ',iI,iA
          call RecPrt('Int1:','(8F10.6)',Int1_2,nOrb(iSym_C)+nDel(iSym_C),nOrb(iSym_D)+nDel(iSym_D))
          write(u6,*) ' *  I,A = ',iI,iA
          call RecPrt('Int2:','(8F10.6)',Int2_2,nOrb(iSym_C)+nDel(iSym_C),nOrb(iSym_D)+nDel(iSym_D))
#         endif
        end if

        ! Construct Tiajb for a specific (i,a)-pair
        if (NonZeroSym(1)) then
          do iJ=1,nJ
            iJJ = iJ+nFro(iSym_C)
            do iB=1,nB
              iBB = iB+nFro(iSym_D)+nOcc(iSym_D)
              xiajb = Int2((iJJ-1)*nTOrb(iSym_D)+iBB)
              xibja = Int1((iJJ-1)*nTOrb(iSym_D)+iBB)
              EDenom = (EOcc(mAdOcc(iSym_A)+iI-1)+EOcc(mAdOcc(iSym_C)+iJ-1)-EVir(mAdVir(iSym_B)+iA-1)-EVir(mAdVir(iSym_D)+iB-1))
              Tiajb = (Two*xiajb-xibja)/EDenom
              Temp1((iJ-1)*nB+iB) = Tiajb
            end do !iSymB
          end do   !iSymJ
#         ifdef _DEBUGPRINT_
          call RecPrt('Tiajb',' ',Temp1,nExt(iSym_D),nOcc(iSym_C))
#         endif
        end if

        if (NonZeroSym(2) .and. (.not. Triangular)) then
          do iJ=1,nJ2
            iJJ = iJ+nFro(iSym_D)
            do iB=1,nB2
              iBB = iB+nFro(iSym_C)+nOcc(iSym_C)
              xiajb = Int2_2((iJJ-1)*nTOrb(iSym_C)+iBB)
              xibja = Int1_2((iJJ-1)*nTOrb(iSym_C)+iBB)
              EDenom = (EOcc(mAdOcc(iSym_A)+iI-1)+EOcc(mAdOcc(iSym_D)+iJ-1)-EVir(mAdVir(iSym_B)+iA-1)-EVir(mAdVir(iSym_C)+iB-1))
              Tiajb = (Two*xiajb-xibja)/EDenom
              Temp1_2((iJ-1)*nB2+iB) = Tiajb
            end do !iSymB
          end do   !iSymJ
#         ifdef _DEBUGPRINT_
          call RecPrt('Tiajb',' ',Temp1,nExt(iSym_D),nOcc(iSym_C))
#         endif
        end if

        ! Do first halftransformation.
        if (NonZeroSym(1)) then
          call dGemm_('N','N',nBas(iSym_D),nOcc(iSym_C),nExt(iSym_D),One,CMO_v(iOffCMO_v(iSym_D)+1),nBas(iSym_D),Temp1, &
                      nExt(iSym_D),Zero,Temp2,nBas(iSym_D))
          call dGemm_('N','T',nBas(iSym_D),nBas(iSym_C),nOcc(iSym_C),One,Temp2,nBas(iSym_D),CMO_o(iOffCMO_o(iSym_C)+1), &
                      nBas(iSym_C),Zero,Temp1,nBas(iSym_D))
        end if

        if (NonZeroSym(2) .and. (.not. Triangular)) then
          call dGemm_('N','N',nBas(iSym_C),nOcc(iSym_D),nExt(iSym_C),One,CMO_v(iOffCMO_v(iSym_C)+1),nBas(iSym_C),Temp1_2, &
                      nExt(iSym_C),Zero,Temp2_2,nBas(iSym_C))
          call dGemm_('N','T',nBas(iSym_C),nBas(iSym_D),nOcc(iSym_D),One,Temp2_2,nBas(iSym_C),CMO_o(iOffCMO_o(iSym_D)+1), &
                      nBas(iSym_D),Zero,Temp1_2,nBas(iSym_C))
        end if
        if (Triangular) then
          ! Do a halfsymmetrization with respect to the two AO-indices if they are the same symmetry.
          do iKap=1,nBas(iSym_D)
            do iLam=1,nBas(iSym_D)
              iLamKap1 = iLam+(iKap-1)*nBas(iSym_D)
              iLamKap2 = iKap+(iLam-1)*nBas(iSym_C)
              Temp2(iLamKap1) = (Temp1(iLamKap1)+Temp1(iLamKap2))/2
            end do
          end do
        else
          ! Do a halfsymmetrization with respect to the two AO-indices if they are different symmetries.
          do iKap=1,nBas(iSym_C)
            do iLam=1,nBas(iSym_D)
              iLamKap1 = iLam+(iKap-1)*nBas(iSym_D)
              iLamKap2 = iKap+(iLam-1)*nBas(iSym_C)
              if (NonZeroSym(1) .and. NonZeroSym(2)) then
                Temp2(iLamKap1) = (Temp1(iLamKap1)+Temp1_2(iLamKap2))/2
              else if (NonZeroSym(1)) then
                Temp2(iLamKap1) = Temp1(iLamKap1)/2
              else if (NonZeroSym(2)) then
                Temp2(iLamKap1) = Temp1_2(iLamKap2)/2
              end if
            end do
          end do
        end if

        ! Place the result in a bin. A bin has fixed Lambda and Kappa and
        ! includes Ti,a,lam,kap as well as i,a-adress
        iIA = iA-1+(iI-1)*nA
        do iKap=1,nBas(iSym_C)
          nLam = nBas(iSym_D)
          if (Triangular) nLam = iKap
          do iLam=1,nLam
            iRec = iLam+(iKap-1)*nBas(iSym_D)
            if (Triangular) then
              iBin = max(iKap,iLam)*(max(iKap,iLam)-1)/2+min(iKap,iLam)
            else
              iBin = iRec
            end if
            !write(u6,*) 'iKap,iLam,iBin=',iKap,iLam,iBin
            iNextX = nint(Bin(1,0,iBin))
            iNextX = iNextX+1
            iBinSize = iNextX*2+2
            if (iBinSize > iBinLength) then
              iLastAdr = iAdrBin
              call dDaFile(LuBin,1,Bin(:,:,iBin),iBinLength,iAdrBin)
              iNextX = 1
              Bin(1,0,iBin) = Zero
              Bin(2,0,iBin) = real(iLastAdr,kind=wp)
            end if
            !write(u6,*) Temp2(iRec),real(iIA,kind=wp),iBin
            Bin(2,iNextX,iBin) = Temp2(iRec)
            Bin(1,iNextX,iBin) = real(iIA,kind=wp)
            Bin(1,0,iBin) = real(iNextX,kind=wp)
          end do
        end do
      end do !iA
    end do   !iI

    !Now there is a second loop over everything to fix the mu-nu-symmetrization

    if (.not. Triangular) then

      do iI=1,nI2
        do iA=1,nA2
          if (NonZeroSym(3)) then
            call Exch(iSym_D,iSym_B,iSym_C,iSym_A,iI+nFro(iSym_B),iA+nFro(iSym_A)+nOcc(iSym_A),Int1,Scr1)
            call Coul(iSym_D,iSym_C,iSym_B,iSym_A,iI+nFro(iSym_B),iA+nFro(iSym_A)+nOcc(iSym_A),Int2,Scr1)
#           ifdef _DEBUGPRINT_
            write(u6,*) ' *  I,A = ',iI,iA
            call RecPrt('Int1_lap2:','(8F10.6)',Int1,nOrb(iSym_D)+nDel(iSym_D),nOrb(iSym_C)+nDel(iSym_C))
            write(u6,*) ' *  I,A = ',iI,iA
            call RecPrt('Int2_lap2:','(8F10.6)',Int2,nOrb(iSym_D)+nDel(iSym_D),nOrb(iSym_C)+nDel(iSym_C))
#           endif
          end if

          if (NonZeroSym(4)) then
            call Exch(iSym_C,iSym_B,iSym_D,iSym_A,iI+nFro(iSym_B),iA+nFro(iSym_A)+nOcc(iSym_A),Int1_2,Scr1)
            call Coul(iSym_C,iSym_D,iSym_B,iSym_A,iI+nFro(iSym_B),iA+nFro(iSym_A)+nOcc(iSym_A),Int2_2,Scr1)
#           ifdef _DEBUGPRINT_
            write(u6,*) ' *  I,A = ',iI,iA
            call RecPrt('Int1_lap2:','(8F10.6)',Int1_2,nOrb(iSym_C)+nDel(iSym_C),nOrb(iSym_D)+nDel(iSym_D))
            write(u6,*) ' *  I,A = ',iI,iA
            call RecPrt('Int2_lap:','(8F10.6)',Int2_2,nOrb(iSym_C)+nDel(iSym_C),nOrb(iSym_D)+nDel(iSym_D))
#           endif
          end if
          ! Construct Tiajb for a specific (i,a)-pair
          if (NonZeroSym(3)) then
            do iJ=1,nJ
              iJJ = iJ+nFro(iSym_C)
              do iB=1,nB
                iBB = iB+nFro(iSym_D)+nOcc(iSym_D)
                xiajb = Int2((iJJ-1)*nTOrb(iSym_D)+iBB)
                xibja = Int1((iJJ-1)*nTOrb(iSym_D)+iBB)
                EDenom = (EOcc(mAdOcc(iSym_B)+iI-1)+EOcc(mAdOcc(iSym_C)+iJ-1)-EVir(mAdVir(iSym_A)+iA-1)-EVir(mAdVir(iSym_D)+iB-1))
                Tiajb = (Two*xiajb-xibja)/EDenom
                Temp1((iJ-1)*nB+iB) = Tiajb
              end do !iSymB
            end do   !iSymJ
          end if

          if (NonZeroSym(4)) then
            do iJ=1,nJ2
              iJJ = iJ+nFro(iSym_D)
              do iB=1,nB2
                iBB = iB+nFro(iSym_C)+nOcc(iSym_C)
                xiajb = Int2_2((iJJ-1)*nTOrb(iSym_C)+iBB)
                xibja = Int1_2((iJJ-1)*nTOrb(iSym_C)+iBB)
                EDenom = (EOcc(mAdOcc(iSym_B)+iI-1)+EOcc(mAdOcc(iSym_D)+iJ-1)-EVir(mAdVir(iSym_A)+iA-1)-EVir(mAdVir(iSym_C)+iB-1))
                Tiajb = (Two*xiajb-xibja)/EDenom

                Temp1_2((iJ-1)*nB2+iB) = Tiajb
              end do !iSymB
            end do   !iSymJ
          end if
          ! Do first halftransformation.
          if (NonZeroSym(3)) then
            call dGemm_('N','N',nBas(iSym_D),nOcc(iSym_C),nExt(iSym_D),One,CMO_v(iOffCMO_v(iSym_D)+1),nBas(iSym_D),Temp1, &
                        nExt(iSym_D),Zero,Temp2,nBas(iSym_D))
            call dGemm_('N','T',nBas(iSym_D),nBas(iSym_C),nOcc(iSym_C),One,Temp2,nBas(iSym_D),CMO_o(iOffCMO_o(iSym_C)+1), &
                        nBas(iSym_C),Zero,Temp1,nBas(iSym_D))
          end if

          if (NonZeroSym(4)) then
            call dGemm_('N','N',nBas(iSym_C),nOcc(iSym_D),nExt(iSym_C),One,CMO_v(iOffCMO_v(iSym_C)+1),nBas(iSym_C),Temp1_2, &
                        nExt(iSym_C),Zero,Temp2_2,nBas(iSym_C))
            call dGemm_('N','T',nBas(iSym_C),nBas(iSym_D),nOcc(iSym_D),One,Temp2_2,nBas(iSym_C),CMO_o(iOffCMO_o(iSym_D)+1), &
                        nBas(iSym_D),Zero,Temp1_2,nBas(iSym_C))
          end if
          ! Do a halfsymmetrization with respect to the two AO-indices if they are the same symmetry.
          do iKap=1,nBas(iSym_C)
            do iLam=1,nBas(iSym_D)
              iLamKap1 = iLam+(iKap-1)*nBas(iSym_D)
              iLamKap2 = iKap+(iLam-1)*nBas(iSym_C)
              if (NonZeroSym(3) .and. NonZeroSym(4)) then
                Temp2(iLamKap1) = (Temp1(iLamKap1)+Temp1_2(iLamKap2))/2
              else if (NonZeroSym(3)) then
                Temp2(iLamKap1) = Temp1(iLamKap1)/2
              else if (NonZeroSym(4)) then
                Temp2(iLamKap1) = Temp1_2(iLamKap2)/2
              end if
            end do
          end do

          iIA = iA-1+(iI-1)*nA2
          do iKap=1,nBas(iSym_C)
            do iLam=1,nBas(iSym_D)
              iRec = iLam+(iKap-1)*nBas(iSym_D)
              iBin = iRec
              !write(u6,*) 'iKap,iLam,iBin=',iKap,iLam,iBin
              iNextX = nint(Bin2(1,0,iBin))
              iNextX = iNextX+1

              iBinSize = iNextX*2+2
              if (iBinSize > iBinLength) then
                call dDaFile(LuBin,1,Bin2(:,:,iBin),iBinLength,iAdrBin)
                iNextX = 1
                Bin2(1,0,iBin) = Zero
                Bin2(2,0,iBin) = real(iAdrBin,kind=wp)
              end if

              !write(u6,*) Temp2(iRec),real(iIA,kind=wp),iBin
              Bin2(2,iNextX,iBin) = Temp2(iRec)
              Bin2(1,iNextX,iBin) = real(iIA,kind=wp)
              Bin2(1,0,iBin) = real(iNextX,kind=wp)
            end do
          end do
        end do !iA
      end do   !iI

    end if

  end if

  ! Read the halftransformed integrals from one bin at the time
  ! and do the two remaining transformations.

  !write(u6,*) 'Do the second half transformation'
  do iBin=1,nBins
    if (NonZeroSym(1) .or. NonZeroSym(2)) then
      do
        iLen = nint(Bin(1,0,iBin))
        do i=1,iLen
          iIA = nint(Bin(1,i,iBin))+1
          Temp1(iIA) = Bin(2,i,iBin)
        end do

        iNextAdr = nint(Bin(2,0,iBin))
        if (iNextAdr == -1) exit
        iAdrRdBin = iNextAdr
        call dDaFile(LuBin,2,Bin(:,:,iBin),iBinLength,iAdrRdBin)
      end do
      ! Do second halftransformation
      call dGemm_('N','N',nBas(iSym_B),nOcc(iSym_A),nExt(iSym_B),One,CMO_v(iOffCMO_v(iSym_B)+1),nBas(iSym_B),Temp1,nExt(iSym_B), &
                  Zero,Temp2,nBas(iSym_B))
      call dGemm_('N','T',nBas(iSym_B),nBas(iSym_A),nOcc(iSym_A),One,Temp2,nBas(iSym_B),CMO_o(iOffCMO_o(iSym_A)+1),nBas(iSym_A), &
                  Zero,Temp1,nBas(iSym_B))
      !call RecPrt('(m,n,i,k)',' ',Temp1,nBas(iSym_B),nBas(iSym_A))
    end if

    if ((NonZeroSym(3) .or. NonZeroSym(4)) .and. (.not. Triangular)) then
      do
        iLen = nint(Bin2(1,0,iBin))
        do i=1,iLen
          iIA = nint(Bin2(1,i,iBin))+1
          Temp1_2(iIA) = Bin2(2,i,iBin)
        end do

        iNextAdr = nint(Bin2(2,0,iBin))
        if (iNextAdr == -1) exit
        iAdrRdBin = iNextAdr
        call dDaFile(LuBin,2,Bin2(:,:,iBin),iBinLength,iAdrRdBin)
      end do
      call dGemm_('N','N',nBas(iSym_A),nOcc(iSym_B),nExt(iSym_A),One,CMO_v(iOffCMO_v(iSym_A)+1),nBas(iSym_A),Temp1_2,nExt(iSym_A), &
                  Zero,Temp2_2,nBas(iSym_A))
      call dGemm_('N','T',nBas(iSym_A),nBas(iSym_B),nOcc(iSym_B),One,Temp2_2,nBas(iSym_A),CMO_o(iOffCMO_o(iSym_B)+1),nBas(iSym_B), &
                  Zero,Temp1_2,nBas(iSym_A))
    end if
    ! Do the second symmetrization. (of the mu and nu-index)
    if (Triangular) then
      do iMu=1,nBas(iSym_A)
        nNu = nBas(iSym_B)
        if (Triangular) nNu = iMu
        do iNu=1,nNu
          iTriMuNu = max(iNu,iMu)*(max(iNu,iMu)-1)/2+min(iNu,iMu)
          iMuNu1 = iNu+(iMu-1)*nBas(iSym_B)
          iMuNu2 = iMu+(iNu-1)*nBas(iSym_A)
          Temp2(iTriMuNu) = (Temp1(iMuNu1)+Temp1(iMuNu2))/2
        end do
      end do
    else
      do iMu=1,nBas(iSym_A)
        do iNu=1,nBas(iSym_B)
          iMuNu1 = iNu+(iMu-1)*nBas(iSym_B)
          iMuNu2 = iMu+(iNu-1)*nBas(iSym_A)
          if ((NonZeroSym(1) .or. NonZeroSym(2)) .and. (NonZeroSym(3) .or. NonZeroSym(4))) then
            Temp2(iMuNu1) = (Temp1(iMuNu1)+Temp1_2(iMuNu2))/2
          else if (NonZeroSym(1) .or. NonZeroSym(2)) then
            Temp2(iMuNu1) = Temp1(iMuNu1)/2
          else if (NonZeroSym(3) .or. NonZeroSym(4)) then
            Temp2(iMuNu1) = Temp1_2(iMuNu2)/2
          end if
        end do
      end do
    end if
    ! Fix the special case with empty iajb:s but "nonempty" (mu nu|kap lam)
    if (LoadZeros) then
      Temp2(:) = Zero
    end if

    ! Store the result on disk.
    if (Triangular) then
      iSize = nBas(iSym_A)*(nBas(iSym_A)+1)/2
#     ifdef _DEBUGPRINT_
      call TriPrt('(Gamma)',' ',Temp2,nBas(iSym_A))
#     endif
    else
      iSize = nBas(iSym_A)*nBas(iSym_B)
#     ifdef _DEBUGPRINT_
      call RecPrt('(Gamma)',' ',Temp2,nBas(iSym_A),nBas(iSym_B))
#     endif
    end if
    call dDaFile(LuGam,1,Temp2,iSize,iAdrGam)
  end do
  LoadZeros = .false.
end do

call DaClos(LuGam)
call DaClos(LuBin)

! Deallocate memory

call mma_deallocate(iTable)
call mma_deallocate(CMO_o)
call mma_deallocate(CMO_v)

call mma_deallocate(Temp1)
call mma_deallocate(Temp2)
call mma_deallocate(Bin)

if (nBlocks /= 1) then
  call mma_deallocate(Temp1_2)
  call mma_deallocate(Temp2_2)
  call mma_deallocate(Bin2)
end if

return

end subroutine Gamma_new
