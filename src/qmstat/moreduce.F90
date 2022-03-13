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

subroutine MoReduce(nBas,MOsToKeep,ipAvRedMO)

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6, r8

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
#include "lenin.fh"
integer(kind=iwp) :: nBas(MxSym), MOsToKeep, ipAvRedMO
integer(kind=iwp) :: i, iAUX, iB, iB1, iB2, icomp, iDav, iDavS, iDin, iDiskUt, iDsq, iM1, iM2, ind, ind1, ind2, ind3, indx, &
                     iNewOcc, iOcc, iopt, iOtD, iOtDt, ipInv, ipTEMP, ipTmoD, ipTreD, ipTreT, irc, iS, iS1, iS2, Ising, iSmLbl, &
                     iSs, iSsq, iSst, iSt, iStrans, iStri, iSx, iTocBig(MxStOT), iTrans, iTransB, iVecs, iVecs2, j, kaunt, & !IFG
                     kaunter, Lu_One, Lu_Scratch, nMtK, nSize
real(kind=wp) :: ChargeNonReduced, ChargeReduced, Det, DiffMax, DiffMegaMax, Dum, Dummy(1), Fac, Sqroot, ThrOcc, TraceFull, &
                 TraceRed, weight
logical(kind=iwp) :: First = .true., LindMOs(MxBas) !IFG
character(len=LenIn8) :: BsLbl(MxBas) !IFG
character(len=50) :: Header
real(kind=wp), parameter :: ReduceWarning = Half
integer(kind=iwp), external :: IsFreeUnit
real(kind=r8), external :: Ddot_

! A word of welcome.

write(u6,*) '     ----- Transform to average natural MO-reduced basis.'

! First we accumulate the different density matrices.

nSize = nBas(1)*(nBas(1)+1)/2
weight = One/real(nState,kind=wp)
call GetMem('DenM','Allo','Real',iDin,nSize)
call GetMem('DenA','Allo','Real',iDav,nSize)
call dcopy_(nSize,[Zero],0,Work(iDav),1)
do iS1=1,nState
  do iS2=1,iS1
    indx = (iS1*(iS1-1)/2+iS2-1)*nSize
    call dcopy_(nSize,Work(iBigT+indx),1,Work(iDin),1)
    if (iS1 /= iS2) cycle
    kaunt = 0
    do iB1=1,nBas(1)
      do iB2=1,iB1
        if (iB1 == iB2) then
          Fac = weight
        else
          Fac = Half*weight
        end if
        Work(iDav+kaunt) = Work(iDav+kaunt)+Work(iDin+kaunt)*Fac
        kaunt = kaunt+1
      end do
    end do
  end do
end do

! Then since we are working in the non-orthogonal AO-basis,
! it is necessary to orthogonalize before we diagonalize
! accumulated density.

call GetMem('Vecs','Allo','Real',iVecs,nBas(1)**2)
call GetMem('Vecs2','Allo','Real',iVecs2,nBas(1)**2)
call GetMem('AuxS','Allo','Real',iAUX,nBas(1)**2)
call GetMem('OvlSs','Allo','Real',iSs,nBas(1)**2)
call GetMem('OvlSs','Allo','Real',iSx,nSize)
call GetMem('OvlSi','Allo','Real',iSst,nBas(1)**2)
call GetMem('OvlSi','Allo','Real',iSt,nSize)
call GetMem('DavSq','Allo','Real',iDavS,nBas(1)**2)
call GetMem('Trans','Allo','Real',iTrans,nBas(1)**2)
call GetMem('Trans','Allo','Real',iTransB,nBas(1)**2)
call GetMem('OrtoAvDen','Allo','Real',iOtD,nBas(1)**2)
call GetMem('OrtoAcDeT','Allo','Real',iOtDt,nSize)
call GetMem('OvlS','Allo','Real',iS,nSize+4)
call GetMem('Occs','Allo','Real',iOcc,nBas(1))
kaunter = 0
do iB1=1,nBas(1)
  do iB2=1,nBas(1)
    Work(iVecs+kaunter) = 0
    if (iB1 == iB2) Work(iVecs+kaunter) = 1
    kaunter = kaunter+1
  end do
end do
! Symmetric orthogonalization, hence get overlap matrix, S.
Lu_One = 92
call OpnOne(irc,0,'ONEINT',Lu_One)
irc = -1
iopt = 0
iSmLbl = 0
icomp = 1
call RdOne(irc,iopt,'Mltpl  0',icomp,Work(iS),iSmLbl)
call Jacob(Work(iS),Work(iVecs),nBas(1),nBas(1))
call dcopy_(nSize,[Zero],0,Work(iSx),1)
call dcopy_(nSize,[Zero],0,Work(iSt),1)
do i=1,nBas(1)
  Sqroot = sqrt(Work(iS+i*(i+1)/2-1))
  Work(iSx+i*(i+1)/2-1) = One/Sqroot
  Work(iSt+i*(i+1)/2-1) = Sqroot
end do
call Square(Work(iSx),Work(iSs),1,nBas(1),nBas(1))
call Square(Work(iSt),Work(iSst),1,nBas(1),nBas(1))
! S^(-1/2)
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Work(iVecs),nBas(1),Work(iSs),nBas(1),Zero,Work(iAUX),nBas(1))
call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),One,Work(iAUX),nBas(1),Work(iVecs),nBas(1),Zero,Work(iTrans),nBas(1))
! S^(1/2)
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Work(iVecs),nBas(1),Work(iSst),nBas(1),Zero,Work(iAUX),nBas(1))
call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),One,Work(iAUX),nBas(1),Work(iVecs),nBas(1),Zero,Work(iTransB),nBas(1))
! The density matrix transforms 'inversely' from the matrix-elements,
! thus let S^(1/2) transform it, not S^(-1/2) which applies to the
! matrix elements.
call Square(Work(iDav),Work(iDavS),1,nBas(1),nBas(1))
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Work(iTransB),nBas(1),Work(iDavS),nBas(1),Zero,Work(iAUX),nBas(1))
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Work(iAUX),nBas(1),Work(iTransB),nBas(1),Zero,Work(iOtD),nBas(1))
kaunter = 0
do iB1=1,nBas(1)
  do iB2=1,nBas(1)
    Work(iVecs2+kaunter) = 0
    if (iB1 == iB2) Work(iVecs2+kaunter) = 1
    kaunter = kaunter+1
  end do
end do
call SqToTri_Q(Work(iOtD),Work(iOtDt),nBas(1))
call Jacob(Work(iOtDt),Work(iVecs2),nBas(1),nBas(1))
! With diagonalized density matrix, collect occupation numbers and
! natural orbital coefficients.
call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Work(iTrans),nBas(1),Work(iVecs2),nBas(1),Zero,Work(iAUX),nBas(1))
kaunt = 0
kaunter = 0
do i=1,nBas(1)
  do j=1,i
    if (i == j) then
      Work(iOcc+kaunt) = Work(iOtDt+kaunter)
      kaunt = kaunt+1
    end if
    kaunter = kaunter+1
  end do
end do
TraceFull = 0
do i=1,nBas(1)
  TraceFull = TraceFull+Work(iOcc+i-1)
end do
if (iPrint >= 10) then
  call Get_cArray('Unique Basis Names',BsLbl,LenIn8*nBas(1))
  write(Header,'(A)') 'All average transition density orbitals'
  ThrOcc = -One
  call Primo(Header,.true.,.false.,ThrOcc,Dum,1,nBas(1),nBas(1),BsLbl,Dummy,Work(iOcc),Work(iAUX),-1)
  write(u6,*)
  write(u6,*) '  Trace = ',TraceFull
end if
! Deallocations.
call GetMem('DenM','Free','Real',iDin,nSize)
call GetMem('DenA','Free','Real',iDav,nSize)
call GetMem('Vecs','Free','Real',iVecs,nBas(1)**2)
call GetMem('Vecs2','Free','Real',iVecs2,nBas(1)**2)
call GetMem('OvlSs','Free','Real',iSs,nBas(1)**2)
call GetMem('OvlSs','Free','Real',iSx,nSize)
call GetMem('OvlSi','Free','Real',iSst,nBas(1)**2)
call GetMem('OvlSi','Free','Real',iSt,nSize)
call GetMem('DavSq','Free','Real',iDavS,nBas(1)**2)
call GetMem('Trans','Free','Real',iTrans,nBas(1)**2)
call GetMem('Trans','Free','Real',iTransB,nBas(1)**2)
call GetMem('OrtoAvDen','Free','Real',iOtD,nBas(1)**2)
call GetMem('OrtoAcDeT','Free','Real',iOtDt,nSize)
call GetMem('OvlS','Free','Real',iS,nSize+4)

! Jetzt far wir mal wieder. Reduce MO-basis according to input criterion.

MOsToKeep = 0
do iB=1,nBas(1)
  if (Work(iOcc+iB-1) >= ThrsRedOcc) then
    LindMOs(iB) = .true.
    MOsToKeep = MOsToKeep+1
  else
    LindMOs(iB) = .false.
  end if
end do
nSize = nBas(1)*MOsToKeep
call GetMem('UncleMoe','Allo','Real',ipAvRedMO,nSize)
call GetMem('NewOccs','Allo','Real',iNewOcc,MOsToKeep)
ind2 = 0
ind3 = 0
! Loop to suck-out the nice MOs.
do iB=1,nBas(1)
  if (LindMOs(iB)) then
    ind1 = nBas(1)*(iB-1)
    call dcopy_(nBas(1),Work(iAUX+ind1),1,Work(ipAvRedMO+ind2),1)
    Work(iNewOcc+ind3) = Work(iOcc+iB-1)
    ind2 = ind2+nBas(1)
    ind3 = ind3+1
  end if
end do
TraceRed = 0
do i=1,MOsToKeep
  TraceRed = TraceRed+Work(iNewOcc+i-1)
end do
! Make a trace check.
if ((TraceFull-TraceRed) >= ReduceWarning) then
  write(u6,*)
  write(u6,*) 'WARNING!  With your occupation threshold, the density matrix trace'
  write(u6,*) 'differs by ',TraceFull-TraceRed,'.'
  write(u6,*) 'You should consider lowering the threshold!'
end if
if (iPrint >= 5) then
  call Get_cArray('Unique Basis Names',BsLbl,LenIn8*nBas(1))
  write(Header,'(A)') 'Reduced average orbitals'
  ThrOcc = -One
  call Primo(Header,.true.,.false.,ThrOcc,Dum,1,nBas(1),[MOsToKeep],BsLbl,Dummy,Work(iNewOcc),Work(ipAvRedMO),-1)
  write(u6,*)
  write(u6,*) '  Trace = ',TraceRed,MOsToKeep
end if

! Time to reduce all individual density matrices in the big TDM to
! the reduced MO-basis. Once more, observe that the density
! transforms contravariantly. But sadly, we need to invert the
! full square MO-matrix before we make reductions.

call GetMem('InverseC','Allo','Real',ipInv,nBas(1)**2)
call MInv(Work(iAUX),Work(ipInv),Ising,Det,nBas(1))

! Now all those transformations and density reductions. To check for
! density losses, the overlaps are read. These partially transformed
! transition density matrix is stored in a scratch file.

DiffMegaMax = Zero
nSize = nBas(1)*(nBas(1)+1)/2
nMtK = MOsToKeep*(MOsToKeep+1)/2
call GetMem('Temporary','Allo','Real',ipTEMP,nBas(1)**2)
call GetMem('MOtrDen','Allo','Real',ipTmoD,nBas(1)**2)
call GetMem('MOreDen','Allo','Real',ipTreD,MOsToKeep**2)
call GetMem('MOreDen','Allo','Real',ipTreT,nMtK)
call GetMem('DenM','Allo','Real',iDin,nSize)
call GetMem('DenMsq','Allo','Real',iDsq,nBas(1)**2)
call GetMem('OvlS','Allo','Real',iS,nSize+4)
call GetMem('Ssquare','Allo','Real',iSsq,nBas(1)**2)
call GetMem('Strans','Allo','Real',iStrans,nBas(1)**2)
call GetMem('Stri','Allo','Real',iStri,nMtK)
Lu_Scratch = 57
Lu_Scratch = IsFreeUnit(Lu_Scratch)
call DaName(Lu_Scratch,'TDMSCR')
irc = -1
iopt = 0
iSmLbl = 0
icomp = 1
call RdOne(irc,iopt,'Mltpl  0',icomp,Work(iS),iSmLbl)
call ClsOne(irc,iopt)
iDiskUt = 0
call iDaFile(Lu_Scratch,1,iTocBig,MxStOT,iDiskUt)
do iS1=1,nState
  do iS2=1,iS1
    ! Collect this particular density matrix.
    indx = (iS1*(iS1-1)/2+iS2-1)*nSize
    call dcopy_(nSize,Work(iBigT+indx),1,Work(iDin),1)
    ! Square it and correct the non-diagonal (recall convention)
    call Dsq(Work(iDin),Work(iDsq),1,nBas(1),nBas(1))
    ! Contravariant transformation of density matrix.
    call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Work(ipInv),nBas(1),Work(iDsq),nBas(1),Zero,Work(ipTEMP),nBas(1))
    call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),One,Work(ipTEMP),nBas(1),Work(ipInv),nBas(1),Zero,Work(ipTmoD),nBas(1))
    ! Covariant transformation of overlap matrix.
    call Square(Work(iS),Work(iSsq),1,nBas(1),nBas(1))
    call Dgemm_('T','N',nBas(1),nBas(1),nBas(1),One,Work(iAUX),nBas(1),Work(iSsq),nBas(1),Zero,Work(ipTEMP),nBas(1))
    call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),One,Work(ipTEMP),nBas(1),Work(iAUX),nBas(1),Zero,Work(iStrans),nBas(1))
    ! How much charge ('overlap') is there in this density element?
    ChargeNonReduced = Ddot_(nBas(1)**2,Work(ipTmoD),1,Work(iStrans),1)
    ! Reduction of density matrix and overlap matrix.
    kaunter = 0
    ind1 = 0
    do iM1=1,nBas(1)
      do iM2=1,nBas(1)
        if (LindMOs(iM1) .and. LindMOs(iM2)) then
          Work(ipTreD+ind1) = Work(ipTmoD+kaunter)
          Work(iSsq+ind1) = Work(iStrans+kaunter)
          ind1 = ind1+1
        end if
        kaunter = kaunter+1
      end do
    end do
    ! Triangularize, jetzt!
    kaunter = 0
    do iM1=1,MOsToKeep
      do iM2=1,MOsToKeep
        if (iM1 /= iM2) Work(ipTreD+kaunter) = 2*Work(ipTreD+kaunter)
        kaunter = kaunter+1
      end do
    end do
    call SqToTri_Q(Work(ipTreD),Work(ipTreT),MOsToKeep)
    call SqToTri_Q(Work(iSsq),Work(iStri),MOsToKeep)
    ! Compute total electronic charge of this reduced density.
    ChargeReduced = Ddot_(nMtK,Work(ipTreT),1,Work(iStri),1)
    ! Renormalize to get right charge ('overlap'); to safeguard
    ! against zero overlaps, make check.
    if ((abs(ChargeNonReduced) <= 1.0e-7_wp) .or. (abs(ChargeReduced) <= 1.0e-7_wp)) then
      Fac = One
    else
      Fac = ChargeNonReduced/ChargeReduced
    end if
    kaunter = 0
    do iM1=1,MOsToKeep
      do iM2=1,iM1
        Work(ipTreT+kaunter) = Work(ipTreT+kaunter)*Fac
        kaunter = kaunter+1
      end do
    end do
    ! If sufficient printlevel, show moment modifications.
    if (iPrint >= 10) then
      call MomentMod(ipTreT,ipTmoD,iAUX,MOsToKeep,nBas(1),LindMOs,iS1,iS2,First,DiffMax)
      if (DiffMax > DiffMegaMax) DiffMegaMax = DiffMax
    end if
    ! Add previous disk address to T-o-C.
    ind = iS1*(iS1+1)/2-iS1+iS2
    iTocBig(ind) = iDiskUt
    ! 'Because I will take a giant dump on you!'
    call dDaFile(Lu_Scratch,1,Work(ipTreT),nMtk,iDiskUt)
  end do
end do
! The real table-of-content
iDiskUt = 0
call iDaFile(Lu_Scratch,1,iTocBig,MxStOT,iDiskUt)
! Deallocations and closing.
call GetMem('NewOccs','Free','Real',iNewOcc,MOsToKeep)
call GetMem('InverseC','Free','Real',ipInv,nBas(1)**2)
call GetMem('Temporary','Free','Real',ipTEMP,nBas(1)**2)
call GetMem('MOtrDen','Free','Real',ipTmoD,nBas(1)**2)
call GetMem('MOreDen','Free','Real',ipTreD,MOsToKeep**2)
call GetMem('MOreDen','Free','Real',ipTreT,nMtK)
call GetMem('DenM','Free','Real',iDin,nSize)
call GetMem('DenMsq','Free','Real',iDsq,nBas(1)**2)
call GetMem('OvlS','Free','Real',iS,nSize+4)
call GetMem('Ssquare','Free','Real',iSsq,nBas(1)**2)
call GetMem('Strans','Free','Real',iStrans,nBas(1)**2)
call GetMem('Stri','Free','Real',iStri,nMtK)
call GetMem('AuxS','Free','Real',iAUX,nBas(1)**2)
call GetMem('Occs','Free','Real',iOcc,nBas(1))
call DaClos(Lu_Scratch)

! Report on the reduction.

write(u6,*)
write(u6,90) 'AO-basis ---> MO-basis reduction complete.'
write(u6,91) 'From ',nBas(1),' functions to ',MosToKeep,'.'
write(u6,90) 'Reduced basis renormalized to have same overlap as non-reduced.'
if (iPrint >= 10) then
  write(u6,92) 'Largest dipole difference is ',DiffMegaMax
end if

return

90 format('        ',A)
91 format('        ',A,I3,A,I3,A)
92 format('        ',A,F10.7)

end subroutine MoReduce
