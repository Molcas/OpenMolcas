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

subroutine MoReduce(nBas,MOsToKeep)

use qmstat_global, only: AvRed, BigT, iPrint, MxSymQ, nState, ThrsRedOcc
use Index_Functions, only: iTri, nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas(MxSymQ)
integer(kind=iwp), intent(out) :: MOsToKeep
#include "Molcas.fh"
integer(kind=iwp) :: i, iB, iB1, iB2, icomp, iDiskUt, iM1, iM2, ind, ind1, ind2, indx, iopt, irc, iS1, iS2, iSmLbl, kaunt, Lu_One, &
                     Lu_Scratch, nMtK, nSize
real(kind=wp) :: ChargeNonReduced, ChargeReduced, Det, DiffMax, DiffMegaMax, Dum, Dummy(1), Sqroot, ThrOcc, TraceFull, TraceRed, &
                 weight
logical(kind=iwp) :: First = .true.
character(len=50) :: Header
character(len=8) :: Label
integer(kind=iwp), allocatable :: iTocBig(:)
real(kind=wp), allocatable :: AUX(:,:), Dav(:), DavS(:,:), Din(:), DsqM(:,:), Inv(:,:), NewOcc(:), Occ(:), OtD(:,:), OtDt(:), &
                              S(:), Ss(:,:), Ssq(:), Sst(:,:), St(:), Strans(:,:), Stri(:), Sx(:), TEMP(:,:), TmoD(:,:), &
                              Trans(:,:), TransB(:,:), TreD(:,:), TreT(:), Vecs(:,:)
logical(kind=iwp), allocatable :: LindMOs(:)
character(len=LenIn8), allocatable :: BsLbl(:)
real(kind=wp), parameter :: ReduceWarning = Half
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: Ddot_

! A word of welcome.

write(u6,*) '     ----- Transform to average natural MO-reduced basis.'

! First we accumulate the different density matrices.

nSize = nTri_Elem(nBas(1))
weight = One/real(nState,kind=wp)
call mma_allocate(Din,nSize,label='DenM')
call mma_allocate(Dav,nSize,label='DenA')
Dav(:) = Zero
do iS1=1,nState
  do iS2=1,iS1
    indx = iTri(iS1,iS2)
    Din(:) = BigT(:,indx)
    if (iS1 /= iS2) cycle
    kaunt = 0
    do iB1=1,nBas(1)
      do iB2=1,iB1-1
        kaunt = kaunt+1
        Dav(kaunt) = Dav(kaunt)+Din(kaunt)*Half*weight
      end do
      kaunt = kaunt+1
      Dav(kaunt) = Dav(kaunt)+Din(kaunt)*weight
    end do
  end do
end do

! Then since we are working in the non-orthogonal AO-basis,
! it is necessary to orthogonalize before we diagonalize
! accumulated density.

call mma_allocate(LindMOs,nBas(1),label='LindMOs')
call mma_allocate(Vecs,nBas(1),nBas(1),label='Vecs')
call mma_allocate(AUX,nBas(1),nBas(1),label='AuxS')
call mma_allocate(Ss,nBas(1),nBas(1),label='OvlSs')
call mma_allocate(Sx,nSize,label='OvlSs')
call mma_allocate(Sst,nBas(1),nBas(1),label='OvlSi')
call mma_allocate(St,nSize,label='OvlSi')
call mma_allocate(DavS,nBas(1),nBas(1),label='DavSq')
call mma_allocate(Trans,nBas(1),nBas(1),label='Trans')
call mma_allocate(TransB,nBas(1),nBas(1),label='Trans')
call mma_allocate(OtD,nBas(1),nBas(1),label='OrtoAvDen')
call mma_allocate(OtDt,nSize,label='OrtoAcDeT')
call mma_allocate(S,nSize,label='OvlS')
call mma_allocate(Occ,nBas(1),label='Occs')
call unitmat(Vecs,nBas(1))
! Symmetric orthogonalization, hence get overlap matrix, S.
Lu_One = 92
iopt = 0
call OpnOne(irc,iopt,'ONEINT',Lu_One)
irc = -1
iopt = ibset(ibset(0,sNoOri),sNoNuc)
iSmLbl = 0
Label = 'Mltpl  0'
icomp = 1
call RdOne(irc,iopt,Label,icomp,S,iSmLbl)
call Jacob(S,Vecs,nBas(1),nBas(1))
Sx(:) = Zero
St(:) = Zero
do i=1,nBas(1)
  Sqroot = sqrt(S(nTri_Elem(i)))
  Sx(nTri_Elem(i)) = One/Sqroot
  St(nTri_Elem(i)) = Sqroot
end do
call Square(Sx,Ss,1,nBas(1),nBas(1))
call Square(St,Sst,1,nBas(1),nBas(1))
! S^(-1/2)
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Vecs,nBas(1),Ss,nBas(1),Zero,AUX,nBas(1))
call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),One,AUX,nBas(1),Vecs,nBas(1),Zero,Trans,nBas(1))
! S^(1/2)
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Vecs,nBas(1),Sst,nBas(1),Zero,AUX,nBas(1))
call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),One,AUX,nBas(1),Vecs,nBas(1),Zero,TransB,nBas(1))
! The density matrix transforms 'inversely' from the matrix-elements,
! thus let S^(1/2) transform it, not S^(-1/2) which applies to the
! matrix elements.
call Square(Dav,DavS,1,nBas(1),nBas(1))
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,TransB,nBas(1),DavS,nBas(1),Zero,AUX,nBas(1))
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,AUX,nBas(1),TransB,nBas(1),Zero,OtD,nBas(1))
call unitmat(Vecs,nBas(1))
call SqToTri_Q(OtD,OtDt,nBas(1))
call Jacob(OtDt,Vecs,nBas(1),nBas(1))
! With diagonalized density matrix, collect occupation numbers and
! natural orbital coefficients.
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Trans,nBas(1),Vecs,nBas(1),Zero,AUX,nBas(1))
do i=1,nBas(1)
  Occ(i) = OtDt(nTri_Elem(i))
end do
TraceFull = Zero
do i=1,nBas(1)
  TraceFull = TraceFull+Occ(i)
end do
if (iPrint >= 5) then
  call mma_allocate(BsLbl,nBas(1),label='BsLbl')
  call Get_cArray('Unique Basis Names',BsLbl,LenIn8*nBas(1))
end if
if (iPrint >= 10) then
  write(Header,'(A)') 'All average transition density orbitals'
  ThrOcc = -One
  call Primo(Header,.true.,.false.,ThrOcc,Dum,1,nBas(1),nBas(1),BsLbl,Dummy,Occ,AUX,-1)
  write(u6,*)
  write(u6,*) '  Trace = ',TraceFull
end if
! Deallocations.
call mma_deallocate(Din)
call mma_deallocate(Dav)
call mma_deallocate(Vecs)
call mma_deallocate(Ss)
call mma_deallocate(Sx)
call mma_deallocate(Sst)
call mma_deallocate(St)
call mma_deallocate(DavS)
call mma_deallocate(Trans)
call mma_deallocate(TransB)
call mma_deallocate(OtD)
call mma_deallocate(OtDt)
call mma_deallocate(S)

! Jetzt far wir mal wieder. Reduce MO-basis according to input criterion.

MOsToKeep = 0
do iB=1,nBas(1)
  if (Occ(iB) >= ThrsRedOcc) then
    LindMOs(iB) = .true.
    MOsToKeep = MOsToKeep+1
  else
    LindMOs(iB) = .false.
  end if
end do
call mma_allocate(AvRed,nBas(1),MOsToKeep,label='UncleMoe')
call mma_allocate(NewOcc,MOsToKeep,label='NewOccs')
ind1 = 0
! Loop to suck-out the nice MOs.
do iB=1,nBas(1)
  if (LindMOs(iB)) then
    ind1 = ind1+1
    AvRed(:,ind1) = AUX(:,iB)
    NewOcc(ind1) = Occ(iB)
  end if
end do
TraceRed = Zero
do i=1,MOsToKeep
  TraceRed = TraceRed+NewOcc(i)
end do
! Make a trace check.
if ((TraceFull-TraceRed) >= ReduceWarning) then
  write(u6,*)
  write(u6,*) 'WARNING!  With your occupation threshold, the density matrix trace'
  write(u6,*) 'differs by ',TraceFull-TraceRed,'.'
  write(u6,*) 'You should consider lowering the threshold!'
end if
if (iPrint >= 5) then
  write(Header,'(A)') 'Reduced average orbitals'
  ThrOcc = -One
  call Primo(Header,.true.,.false.,ThrOcc,Dum,1,nBas(1),[MOsToKeep],BsLbl,Dummy,NewOcc,AvRed,-1)
  write(u6,*)
  write(u6,*) '  Trace = ',TraceRed,MOsToKeep
  call mma_deallocate(BsLbl)
end if

! Time to reduce all individual density matrices in the big TDM to
! the reduced MO-basis. Once more, observe that the density
! transforms contravariantly. But sadly, we need to invert the
! full square MO-matrix before we make reductions.

call mma_allocate(Inv,nBas(1),nBas(1),label='InverseC')
call MInv(AUX,Inv,Det,nBas(1))

! Now all those transformations and density reductions. To check for
! density losses, the overlaps are read. These partially transformed
! transition density matrix is stored in a scratch file.

DiffMegaMax = Zero
nSize = nTri_Elem(nBas(1))
nMtK = nTri_Elem(MOsToKeep)
call mma_allocate(TEMP,nBas(1),nBas(1),label='Temporary')
call mma_allocate(TmoD,nBas(1),nBas(1),label='MOtrDen')
call mma_allocate(TreD,MOsToKeep,MOsToKeep,label='MOreDen')
call mma_allocate(TreT,nMtK,label='MOreDen')
call mma_allocate(Din,nSize,label='DenM')
call mma_allocate(DsqM,nBas(1),nBas(1),label='DenMsq')
call mma_allocate(S,nSize,label='OvlS')
call mma_allocate(Ssq,nBas(1)**2,label='Ssquare')
call mma_allocate(Strans,nBas(1),nBas(1),label='Strans')
call mma_allocate(Stri,nMtK,label='Stri')
Lu_Scratch = IsFreeUnit(57)
call DaName(Lu_Scratch,'TDMSCR')
irc = -1
iopt = ibset(ibset(0,sNoOri),sNoNuc)
iSmLbl = 0
icomp = 1
call RdOne(irc,iopt,Label,icomp,S,iSmLbl)
call ClsOne(irc,iopt)
iDiskUt = 0
call mma_allocate(iTocBig,nTri_Elem(nState),label='iTocBig')
call iDaFile(Lu_Scratch,1,iTocBig,nTri_Elem(nState),iDiskUt)
do iS1=1,nState
  do iS2=1,iS1
    ! Collect this particular density matrix.
    indx = iTri(iS1,iS2)
    Din(:) = BigT(:,indx)
    ! Square it and correct the non-diagonal (recall convention)
    call Dsq(Din,DsqM,1,nBas(1),nBas(1))
    ! Contravariant transformation of density matrix.
    call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Inv,nBas(1),DsqM,nBas(1),Zero,TEMP,nBas(1))
    call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),One,TEMP,nBas(1),Inv,nBas(1),Zero,TmoD,nBas(1))
    ! Covariant transformation of overlap matrix.
    call Square(S,Ssq,1,nBas(1),nBas(1))
    call Dgemm_('T','N',nBas(1),nBas(1),nBas(1),One,AUX,nBas(1),Ssq,nBas(1),Zero,TEMP,nBas(1))
    call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,TEMP,nBas(1),AUX,nBas(1),Zero,Strans,nBas(1))
    ! How much charge ('overlap') is there in this density element?
    ChargeNonReduced = Ddot_(nBas(1)**2,TmoD,1,Strans,1)
    ! Reduction of density matrix and overlap matrix.
    ind1 = 0
    do iM1=1,nBas(1)
      if (.not. LindMOs(iM1)) cycle
      ind1 = ind1+1
      ind2 = 0
      do iM2=1,nBas(1)
        if (.not. LindMOs(iM2)) cycle
        ind2 = ind2+1
        TreD(ind2,ind1) = TmoD(iM2,iM1)
        Ssq(ind2+(ind1-1)*MOsToKeep) = Strans(iM2,iM1)
      end do
    end do
    ! Triangularize, jetzt!
    do iM1=1,MOsToKeep
      do iM2=1,MOsToKeep
        if (iM1 /= iM2) TreD(iM2,iM1) = Two*TreD(iM2,iM1)
      end do
    end do
    call SqToTri_Q(TreD,TreT,MOsToKeep)
    call SqToTri_Q(Ssq,Stri,MOsToKeep)
    ! Compute total electronic charge of this reduced density.
    ChargeReduced = Ddot_(nMtK,TreT,1,Stri,1)
    ! Renormalize to get right charge ('overlap'); to safeguard
    ! against zero overlaps, make check.
    if ((abs(ChargeNonReduced) > 1.0e-7_wp) .and. (abs(ChargeReduced) > 1.0e-7_wp)) then
      TreT(:) = TreT*ChargeNonReduced/ChargeReduced
    end if
    ! If sufficient printlevel, show moment modifications.
    if (iPrint >= 10) then
      call MomentMod(TreT,TmoD,AUX,MOsToKeep,nBas(1),LindMOs,iS1,iS2,First,DiffMax)
      if (DiffMax > DiffMegaMax) DiffMegaMax = DiffMax
    end if
    ! Add previous disk address to T-o-C.
    ind = iTri(iS1,iS2)
    iTocBig(ind) = iDiskUt
    ! 'Because I will take a giant dump on you!'
    call dDaFile(Lu_Scratch,1,TreT,nMtk,iDiskUt)
  end do
end do
! The real table-of-content
iDiskUt = 0
call iDaFile(Lu_Scratch,1,iTocBig,nTri_Elem(nState),iDiskUt)
! Deallocations and closing.
call mma_deallocate(iTocBig)
call mma_deallocate(LindMOs)
call mma_deallocate(AUX)
call mma_deallocate(NewOcc)
call mma_deallocate(Inv)
call mma_deallocate(TEMP)
call mma_deallocate(TmoD)
call mma_deallocate(TreD)
call mma_deallocate(TreT)
call mma_deallocate(Din)
call mma_deallocate(DsqM)
call mma_deallocate(S)
call mma_deallocate(Ssq)
call mma_deallocate(Strans)
call mma_deallocate(Stri)
call mma_deallocate(Occ)
call DaClos(Lu_Scratch)

! Report on the reduction.

write(u6,*)
write(u6,90) 'AO-basis ---> MO-basis reduction complete.'
write(u6,91) 'From ',nBas(1),' functions to ',MosToKeep,'.'
write(u6,90) 'Reduced basis renormalized to have same overlap as non-reduced.'
if (iPrint >= 10) write(u6,92) 'Largest dipole difference is ',DiffMegaMax

return

90 format('        ',A)
91 format('        ',A,I3,A,I3,A)
92 format('        ',A,F10.7)

end subroutine MoReduce
