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

subroutine RdJobIph_td(CIVec)
!***********************************************************************
!                                                                      *
!     Read the contents of the JOBIPH file.                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, FnJob, G1t, G2sq, G2t, LuJob, nA, nNA
use input_mclr, only: ERASSCF, Headerjp, iPT2, iRoot, iSpin, iTOC, iTocIph, lRoots, nActEl, nAsh, nBas, nCOnf, nDel, nElec3, nFro, &
                      nHole1, nIsh, nOrb, nRoots, nRS1, nRS2, nRS3, nSym, ntAsh, ntASqr, ntATri, ntBas, ntBSqr, ntBTri, ntIsh, &
                      ntISqr, ntITri, State_Sym, TitleJP, Weight
use Molcas, only: LenIn, MxOrb, MxRoot, MxSym
use RASDim, only: MxIter, MxTit
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), allocatable, intent(out) :: CIVec(:,:)
integer(kind=iwp) :: i, iB, iDij, iDij2, iDisk, iDkl, iDkl2, iIJKL, iIJKL2, iS, iSym, Iter, jB, jS, kB, kRoots, kS, lB, Length, &
                     lS, nAct, nAct2, nAct4, nG1, nG2
real(kind=wp) :: Fact, Fact2, Factij, Factkl, PotNuc0, rdum(1), Temp
real(kind=wp), allocatable :: G2tta(:), G2tts(:), Tmp2(:)
character, allocatable :: TempTxt(:)

!----------------------------------------------------------------------*
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
call WR_RASSCF_Info(LuJob,2,iDisk,nActEl,iSpin,nSym,State_sym,nFro,nIsh,nAsh,nDel,nBas,MxSym,TempTxt,(LenIn+8)*mxorb,nConf, &
                    HeaderJP,144,TitleJP,4*18*mxTit,PotNuc0,lRoots,nRoots,iRoot,mxRoot,nRs1,nRs2,nRs3,nHole1,nElec3,iPt2,Weight)
call mma_deallocate(TempTxt)
!----------------------------------------------------------------------*
!     Overwrite the variable lroots if approriate                      *
!----------------------------------------------------------------------*
if (kRoots /= -1) then
  if (iPt2 /= 0) then
    write(u6,*)
    write(u6,*) ' *** Error in subroutine RDJOBIPH_TD ***'
    write(u6,*) ' Pt2 /= 0'
    write(u6,*)
  else if (kRoots > lRoots) then
    write(u6,*)
    write(u6,*) ' *** Error in subroutine RDJOBIPH_TD ***'
    write(u6,*) ' kRoots > lRoots'
    write(u6,*)
  end if
  lRoots = kRoots
  nRoots = 1
else if (nRoots /= 1) then
  write(u6,*)
  write(u6,*) ' *** Error in subroutine RDJOBIPH_TD ***'
  write(u6,*) ' nRoots /= 1'
  write(u6,*)
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
!----------------------------------------------------------------------*
!     Load the orbitals used in the last macro iteration               *
!----------------------------------------------------------------------*

call mma_allocate(CMO,Length,Label='CMO')
call Get_dArray_chk('Last orbitals',CMO,Length)
!iDisk = iToc(9)
!if (IPT2 == 0) iDisk = iToc(2)
!call dDaFile(LuJob,2,CMO,ntBsqr,iDisk)
!jpCMO = 1
!do iSym=1,nSym
!  write(Line,'(A,i2.2)') 'MO coefficients, iSym = ',iSym
!  call RecPrt(Line,' ',CMO(jpCMO),nBas(iSym),nBas(iSym))
!  jpCMO = jpCMO+nBas(iSym)*nBas(iSym)
!end do
!----------------------------------------------------------------------*
!     Load the CI vector for the root lRoots                           *
!----------------------------------------------------------------------*
call mma_allocate(CIVec,nConf,1,Label='CIVec')
iDisk = iToc(4)
do i=1,lroots-1
  call dDaFile(LuJob,0,CIVec,nConf,iDisk)
end do
call dDaFile(LuJob,2,CIVec,nConf,iDisk)
if (.false.) call DVcPrt('CI coefficients',' ',CIVec,nConf)
!----------------------------------------------------------------------*
!     Load state energy                                                *
!----------------------------------------------------------------------*
call mma_allocate(Tmp2,mxRoot*mxIter,Label='Tmp2')
iDisk = iToc(6)
call dDaFile(LuJob,2,Tmp2,mxRoot*mxIter,iDisk)
ERASSCF(1) = Zero
do iter=0,mxIter-1
  Temp = Tmp2(iter*mxRoot+lRoots)
  if (Temp /= Zero) ERASSCF(1) = Temp
end do
call mma_deallocate(Tmp2)
!if (debug) write(u6,*) ' RASSCF energy =',ERASSCF(1)

nAct = 0
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
!----------------------------------
! One electron dens - triang stor.
!----------------------------------
nG1 = nTri_Elem(nAct)
call mma_allocate(G1t,nG1,Label='G1t')
G1t(:) = Zero

!---------------------------------------
! Triangular part of two electron dens,
! symmetric part
!---------------------------------------
nG2 = nTri_Elem(nG1)

call mma_allocate(G2sq,nAct**4,Label='G2sq')
call mma_allocate(G2t,nG2,Label='G2t')
call mma_allocate(G2tts,nG2,Label='G2tts')
call mma_allocate(G2tta,nG2,Label='G2tta')
iDisk = iToc(3)
do i=1,lroots-1
  call dDaFile(LuJob,0,rdum,nG1,iDisk)
  call dDaFile(LuJob,0,rdum,nG1,iDisk)
  call dDaFile(LuJob,0,rdum,nG2,iDisk)
  call dDaFile(LuJob,0,rdum,nG2,iDisk)
end do
call dDaFile(LuJob,2,G1t,nG1,iDisk)
call dDaFile(LuJob,0,rdum,nG1,iDisk)
call dDaFile(LuJob,2,G2tts,nG2,iDisk)
call dDaFile(LuJob,2,G2tta,nG2,iDisk)

do iB=1,nAct
  do jB=1,iB
    iDij = iTri(ib,jB)
    do kB=1,ib
      do lB=1,kB
        iDkl = iTri(kB,lB)
        fact = One
        if ((iDij >= iDkl) .and. (kB == lB)) fact = Two
        if ((iDij < iDkl) .and. (iB == jB)) fact = Two
        iijkl = iTri(iDij,iDkl)
        G2t(iijkl) = Fact*G2tts(iijkl)
      end do
    end do
  end do
end do

!----------------------------------------------
! Rectangular part of the two el dens,
! symmetric and asymmetric contributions added
! ipG2sq = ipG2tts + ipG2tta
!----------------------------------------------

do iB=1,nAct
  do jB=1,nact
    Factij = One
    if (ib > jb) Factij = -One
    iDij = iTri(ib,jB)
    iDij2 = ib+(jb-1)*NACT
    do kB=1,nact
      do lB=1,nact
        Factkl = One
        if (kb > lb) Factkl = -One
        iDkl = iTri(kB,lB)
        iDkl2 = kb+(lb-1)*NACT
        fact = One
        Fact2 = Factij*Factkl
        if ((iDij >= iDkl) .and. (kB == lB)) fact = Two
        if ((iDij < iDkl) .and. (iB == jB)) fact = Two
        iijkl = iTri(iDij,iDkl)
        iijkl2 = iDij2+nact**2*(iDkl2-1)

        G2sq(iijkl2) = Fact*(G2tts(iijkl)+G2tta(iijkl)*Fact2)
      end do
    end do
  end do
end do
call mma_deallocate(G2tts)
call mma_deallocate(G2tta)
!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*

end subroutine RdJobIph_td
