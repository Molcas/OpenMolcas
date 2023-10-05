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
! Copyright (C) 2010-2012, Jonas Bostrom                               *
!               2013, Roland Lindh                                     *
!***********************************************************************

subroutine ChoMP2g_GradSetup(irc,CMO)
!***********************************************************************
!                                                                      *
!     Author: Jonas Bostrom (2010-2012)                                *
!             Roland Lindh, Implementation of symmetry (some while in  *
!                           New Orleans during the 245th National ACS  *
!                           meeting, 7-11 April 2013).                 *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use Cholesky, only: nBas, nSym, NumCho
use ChoMP2, only: iAdrOff, lUnit_F, nFro, nMoMo, nMP2Vec, nOcc, nOrb, nVir
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Eight, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp) :: i, iAdr, iAdr1, iAdr2, iAdrB(8), iAdrR(8), iB3, iClos, iCMOo, iComp, iiJ, iiK, iJ, iK, indx, iOff, iOff0, &
                     iOff1, iOff2, iOff3, iOff_A, iOff_A2, iOff_JL, iOff_Ria, iOff_Rin, iOff_Rmn, iOff_WJL(8), iOffA, iOffB, &
                     iOffB1, iOffB2, iOffB3(8,8), iOffBB, iOffCInv(8), iOffCMOo(8), iOffD(8), iOffFF, iOffFO, iOffFV, iOffL, &
                     iOffLRb(8,8), iOffLRo(8,8), iOffOF, iOffOO, iOffOV, iOffR, iOffVF, iOffVO, iOffVV, iOffZ, iOpt, iRd, iSeed, &
                     iSym, iSymlbl, iTypL, iTypR, iVecFF, iVecFO, iVecFV, iVecOF, iVecOO, iVecOV, iVecVF, iVecVO, iVecVV, iWr, j, &
                     jAdr, jSym, k, k_o, k_v, kl, kl_s, kSym, l, lA1, lA2, lB1kl, lB2kl, lB3jl, lB3kl, lB3kl_s, lCab, lCai, lCaK, &
                     lCia, lCij, lCiK, lCJK, lCKa, lCKi, lCmn, lCpn, lCpq, lk, lRecDens, lRia, lRin, lRmn, lSym, lTot, lTriDens, &
                     LuA(2), LuB(2), LuBTmp, lUkn, LuLVec, LuRVec, LuXVec, LuYVec, lVkn, lWJL, lWmjKJ, maxvalue, mSym, mSym2, &
                     nB3(8), nBas2, nBB, nJ, nK, nLRb(8), nLRo(8), nOccAll(8), nOccBas, nOO, nOrbBas, nVec, nVirBas
real(kind=wp) :: Fac, Fact
character(len=8) :: label
character(len=6) :: fname2
character(len=5) :: fname
real(kind=wp), allocatable :: A1(:), A2(:), B1kl(:), B2kl(:), B3jl(:), B3kl(:), B3kl_s(:), Cab(:), Cai(:), CaK(:), Cia(:), Cij(:), &
                              CiK(:), CJK(:), CKa(:), CKi(:), Cmn(:), CMO_Inv(:), CMO_o(:), CMO_v(:), Cpn(:), Cpq(:), &
                              MP2Density(:), MP2TDensity(:), MP2TTotDensity(:), Ria(:), Rin(:), Rmn(:), SCFDensity(:), &
                              SCFTDensity(:), SMat(:), STMat(:), Ukn(:), Vkn(:), WJL(:), WmjKJ(:)
character(len=*), parameter :: SecNam = 'ChoMP2g_GradSetup'
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
iWr = 1
iRd = 2
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute some sizes of some arrays

nBas2 = 0
lTriDens = 0
lRecDens = 0
nOrbBas = 0
nOccBas = 0
nVirBas = 0
iAdr = 1
jAdr = 1
iOff_JL = 0
iCMOo = 0
do iSym=1,nSym
  iOffCMOo(iSym) = iCMOo
  iCMOo = iCMOo+nBas(iSym)*nOcc(iSym)
  iOff_WJL(iSym) = iOff_JL
  iOff_JL = iOff_JL+nMP2Vec(iSym)*NumCho(iSym)
  iAdrR(iSym) = iAdr
  iAdrB(iSym) = jAdr
  iOffD(iSym) = nBas2
  nBas2 = nBas2+nBas(iSym)*nBas(iSym)
  lTriDens = lTriDens+nTri_Elem(nBas(iSym))
  lRecDens = lRecDens+nBas(iSym)*nBas(iSym)
  nOccAll(iSym) = nFro(iSym)+nOcc(iSym)
  iOffCInv(iSym) = nOrbBas
  nOrbBas = nOrbBas+nOrb(iSym)*nBas(iSym)
  nOccBas = nOccBas+nOcc(iSym)*nBas(iSym)
  nVirBas = nVirBas+nVir(iSym)*nBas(iSym)
  nOO = 0
  nBB = 0
  iB3 = 0
  do jSym=1,nSym
    kSym = ieor(iSym-1,jSym-1)+1

    iOffLRo(iSym,jSym) = nOO
    iOffLRb(iSym,jSym) = nBB
    nOO = nOO+nOrb(jSym)*nOrb(kSym)
    nBB = nBB+nBas(jSym)*nBas(kSym)
    iOffB3(iSym,jSym) = iB3
    iB3 = iB3+nBas(jSym)*nOcc(kSym)
  end do
  nB3(iSym) = iB3
  nLRo(iSym) = nOO
  nLRb(iSym) = nBB
  iAdr = iAdr+nBB*nMP2Vec(iSym)
  jAdr = jAdr+nBB*NumCho(iSym)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(CMO_inv,nOrbBas,Label='CMO_Inv')
CMO_inv(:) = Zero
call mma_allocate(CMO_o,nOccBas,Label='CMO_o')
call mma_allocate(CMO_v,nVirBas,Label='CMO_v')
!                                                                      *
!***********************************************************************
!                                                                      *
! Get memory for density matricies and overlap matrix

call mma_allocate(MP2Density,lRecDens,Label='MP2Density')
call mma_allocate(SCFDensity,lRecDens,Label='SCFDensity')
call mma_allocate(SMat,lRecDens,Label='SMat')

call mma_allocate(MP2TTotDensity,lTriDens,Label='MP2TTotDensity')
call mma_allocate(MP2TDensity,lTriDens,Label='MP2TDensity')
call mma_allocate(SCFTDensity,lTriDens,Label='SCFTDensity')
call mma_allocate(STMat,lTriDens,Label='STMat')
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the Variational MP2 Density and the HF density

call Get_D1ao_Var(MP2TTotDensity,lTriDens)
call Get_dArray_chk('D1ao',SCFTDensity,lTriDens)
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the overlap matrix

iSymlbl = 1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'Mltpl  0'
iComp = 1
call RdOne(irc,iOpt,Label,iComp,STmat,iSymlbl)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the differential MP2 density
MP2TDensity(:) = Two*(MP2TTotDensity(:)-SCFTDensity(:))
!                                                                      *
!***********************************************************************
!                                                                      *
! Modified for triangular storage to quadratic storage. Some of the
! triangular storage is folded, i.e off diagonal element are
! multiplied with 2.

indx = 0
iOff = 0
do iSym=1,nSym
  do i=1,nBas(iSym)
    do j=1,i-1
      indx = indx+1
      MP2Density(j+nBas(iSym)*(i-1)+iOff) = MP2TDensity(indx)*Half
      MP2Density(i+nBas(iSym)*(j-1)+iOff) = MP2TDensity(indx)*Half
      SCFDensity(j+nBas(iSym)*(i-1)+iOff) = SCFTDensity(indx)*Half
      SCFDensity(i+nBas(iSym)*(j-1)+iOff) = SCFTDensity(indx)*Half
      Smat(j+nBas(iSym)*(i-1)+iOff) = STmat(indx)
      Smat(i+nBas(iSym)*(j-1)+iOff) = STmat(indx)

    end do
    indx = indx+1
    MP2Density(i+nBas(iSym)*(i-1)+iOff) = MP2TDensity(indx)
    SCFDensity(i+nBas(iSym)*(i-1)+iOff) = SCFTDensity(indx)
    Smat(i+nBas(iSym)*(i-1)+iOff) = STmat(indx)

  end do
  iOff = iOff+nBas(iSym)**2
end do

! Deallocate temporary memory

call mma_deallocate(MP2TTotDensity)
call mma_deallocate(MP2TDensity)
call mma_deallocate(SCFTDensity)
call mma_deallocate(STMat)
!                                                                      *
!***********************************************************************
!                                                                      *
! Calculate inverse CMO-matrix as C^-1 = C^T*S
! --------------------------------------------

iOff1 = 1
iOff2 = 1
do iSym=1,nSym
  call dGemm_('T','N',nOrb(iSym),nBas(iSym),nBas(iSym),One,CMO(iOff2),nBas(iSym),Smat(iOff1),nBas(iSym),Zero,CMO_inv(iOff2), &
              nOrb(iSym))

  iOff1 = iOff1+nBas(iSym)**2
  iOff2 = iOff2+nOrb(iSym)*nBas(iSym)
end do

call mma_deallocate(SMat)
!                                                                      *
!***********************************************************************
!                                                                      *
! Calculate occupied and virtual inverse CMO-matrices
! ---------------------------------------------------

k_v = 1
k_o = 1
iOff2 = 0
do iSym=1,nSym
  do i=1,nOrb(iSym)
    do j=1,nBas(iSym)
      if (i > nFro(iSym) .and. (i <= nOccAll(iSym))) then
        CMO_o(k_o) = CMO((i-1)*nBas(iSym)+j+iOff2)
        k_o = k_o+1
      else if (i > nOccAll(iSym)) then
        CMO_v(k_v) = CMO((i-1)*nBas(iSym)+j+iOff2)
        k_v = k_v+1
      end if
    end do
  end do
  iOff2 = iOff2+nOrb(iSym)*nBas(iSym)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Open Cholesky vector files.
! ---------------------------

iTypL = 1
iTypR = 2

do iSym=1,nSym
  call ChoMP2_OpenF(1,iTypL,iSym)
  call ChoMP2_OpenF(1,iTypR,iSym)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Open some temporary files

iSeed = 7
LuXVec = IsFreeUnit(iSeed)
FName = 'TMPV1'
call DaName_MF_WA(LuXVec,Fname)

iSeed = 7
LuYVec = isFreeUnit(iSeed)
FName = 'TMPV2'
call DaName_MF_WA(LuYVec,Fname)

iSeed = 7
LuRVec = isFreeUnit(iSeed)
FName = 'TMPV3'
call DaName_MF_WA(LuRVec,Fname)

LuLVec = isFreeUnit(iSeed)
FName = 'TMPV4'
call DaName_MF_WA(LuLVec,Fname)

iSeed = 7
LuBTmp = isFreeUnit(iSeed)
FName = 'TMPV6'
call DaName_MF_WA(LuBTmp,Fname)
!                                                                      *
!***********************************************************************
!                                                                      *
! Open files for the terms in Eqs. 35, 36, 40, and 41.
! i=1 Coulombic and i=2 Exchange contributions.
! This can be merged into singel 2-center and 3-center
! contributions

do i=1,2

  iSeed = 7
  LuA(i) = isFreeUnit(iSeed)
  write(Fname2,'(A5,I1)') 'AMP2V',i
  call DaName_MF_WA(LuA(i),Fname2)

  iSeed = 7
  LuB(i) = isFreeUnit(iSeed)
  write(Fname2,'(A5,I1)') 'BMP2V',i
  call DaName_MF_WA(LuB(i),Fname2)

end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Max number of vectors in one batch. This is not foolproof, if you
! run the code with too little ram you might get problems
! ------------------------------------------------------------------
!
! This has to be recoded to be dynamic. To be done after the
! symmetry has been implemented in the code.

maxvalue = 10
!                                                                      *
!***********************************************************************
!                                                                      *
nVec = 0
do iSym=1,nSym
  nVec = min(max(NumCho(iSym),nMP2Vec(iSym),nVec),maxvalue)
  if (nVec < 1) call SysAbendMsg(SecNam,'Insufficient memory','[1]')
end do
iSym = 1  ! Here we are!

! Allocate memory for L-vectors
! -----------------------------

iVecFF = 1
iVecOF = 2
iVecVF = 3
iVecFO = 4
iVecOO = 5
iVecVO = 6
iVecFV = 7
iVecOV = 8
iVecVV = 9

lCJK = nMoMo(iSym,iVecFF)*nVec
call mma_allocate(CJK,lCJK,Label='CJK')

lCKi = nMoMo(iSym,iVecOF)*nVec
call mma_allocate(CKi,lCKi,Label='CKi')

lCKa = nMoMo(iSym,iVecVF)*nVec
call mma_allocate(CKa,lCKa,Label='CKa')

lCiK = nMoMo(iSYm,iVecFO)*nVec
call mma_allocate(CiK,lCiK,Label='CiK')

lCij = nMoMo(iSym,iVecOO)*nVec
call mma_allocate(Cij,lCij,Label='Cij')

lCia = nMoMo(iSym,iVecVO)*nVec
call mma_allocate(Cia,lCia,Label='Cia')

lCaK = nMoMo(iSym,iVecFV)*nVec
call mma_allocate(CaK,lCaK,Label='CaK')

lCai = nMoMo(iSym,iVecOV)*nVec
call mma_allocate(Cai,lCai,Label='Cai')

lCab = nMoMo(iSym,iVecVV)*nVec
call mma_allocate(Cab,lCab,Label='Cab')

lCpq = nOrb(iSym)*nOrb(iSym)*nVec
call mma_allocate(Cpq,lCpq,Label='Cpq')

! The Cholesky vectors of the MP2 amplitudes,
! see Eq. 26.

lRia = nOcc(iSym)*nVir(iSym)*nVec
call mma_allocate(Ria,lRia,Label='Ria')

lCpn = nOrb(iSym)*nBas(iSym)*nVec
call mma_allocate(Cpn,lCpn,Label='Cpn')

lRin = nOcc(iSym)*nBas(iSym)*nVec
call mma_allocate(Rin,lRin,Label='Rin')

lCmn = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(Cmn,lCmn,Label='Cmn')

lRmn = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(Rmn,lRmn,Label='Rmn')

lUkn = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(Ukn,lUkn,Label='Ukn')

lVkn = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(Vkn,lVkn,Label='Vkn')

lWJL = NumCho(iSym)*nMP2Vec(iSym)
call mma_allocate(WJL,lWJL,Label='WJL')
WJL(:) = Zero

lWmjKJ = nVec*nVec*nOcc(iSym)*nBas(iSym)
call mma_allocate(WmjKJ,lWmjKJ,Label='WmjKJ')

lB3jl = nBas(iSym)*nOcc(iSym)*nVec
call mma_allocate(B3jl,lB3jl,Label='B3jl')

lB3kl = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(B3kl,lB3kl,Label='B3kl')
!                                                                      *
!***********************************************************************
!                                                                      *
! Start loop over batches of C_pq^J vectors

do iSym=1,nSym
  do iiJ=1,NumCho(iSym),nVec
    nJ = min(nVec,NumCho(iSym)-(iiJ-1))

    ! Read C_pq^J-vectors from disk
    ! -----------------------------

    lTot = nMoMo(iSym,iVecFF)*nJ
    iAdr = 1+nMoMo(iSym,iVecFF)*(iiJ-1)+iAdrOff(iSym,iVecFF)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,CJK,lTot,iAdr)

    lTot = nMoMo(iSym,iVecOF)*nJ
    iAdr = 1+nMoMo(iSym,iVecOF)*(iiJ-1)+iAdrOff(iSym,iVecOF)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,CKi,lTot,iAdr)

    lTot = nMoMo(iSym,iVecVF)*nJ
    iAdr = 1+nMoMo(iSym,iVecVF)*(iiJ-1)+iAdrOff(iSym,iVecVF)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,CKa,lTot,iAdr)

    lTot = nMoMo(iSym,iVecFO)*nJ
    iAdr = 1+nMoMo(iSym,iVecFO)*(iiJ-1)+iAdrOff(iSym,iVecFO)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,CiK,lTot,iAdr)

    lTot = nMoMo(iSym,iVecOO)*nJ
    iAdr = 1+nMoMo(iSym,iVecOO)*(iiJ-1)+iAdrOff(iSym,iVecOO)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,Cij,lTot,iAdr)

    lTot = nMoMo(iSym,iVecVO)*nJ
    iAdr = 1+nMoMo(iSym,iVecVO)*(iiJ-1)+iAdrOff(iSym,iVecVO)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,Cia,lTot,iAdr)

    lTot = nMoMo(iSym,iVecFV)*nJ
    iAdr = 1+nMoMo(iSym,iVecFV)*(iiJ-1)+iAdrOff(iSym,iVecFV)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,CaK,lTot,iAdr)

    lTot = nMoMo(iSym,iVecOV)*nJ
    iAdr = 1+nMoMo(iSym,iVecOV)*(iiJ-1)+iAdrOff(iSym,iVecOV)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,Cai,lTot,iAdr)

    lTot = nMoMo(iSym,iVecVV)*nJ
    iAdr = 1+nMoMo(iSym,iVecVV)*(iiJ-1)+iAdrOff(iSym,iVecVV)
    call dDaFile(lUnit_F(iSym,iTypL),iRd,Cab,lTot,iAdr)

    ! Put the C_pq^J-vectors together to one large vector
    ! ------------------------------------------------

    iOff3 = 0
    do iJ=1,nJ

      iOffFF = 0
      iOffOF = 0
      iOffFO = 0
      iOffOO = 0
      iOffVF = 0
      iOffFV = 0
      iOffVO = 0
      iOffOV = 0
      iOffVV = 0
      iOffBB = 0
      do jSym=1,nSym
        kSym = ieor(iSym-1,jSym-1)+1
        do i=1,nFro(kSym)
          iOff1 = iOffFF+nMoMo(iSym,iVecFF)*(iJ-1)+(i-1)*nFro(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nFro(jSym)) = CJK(iOff1+1:iOff1+nFro(jSym))
          iOff3 = iOff3+nFro(jSym)
          iOff1 = iOffOF+nMoMo(iSym,iVecOF)*(iJ-1)+(i-1)*nOcc(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nOcc(jSym)) = Cki(iOff1+1:iOff1+nOcc(jSym))
          iOff3 = iOff3+nOcc(jSym)
          iOff1 = iOffVF+nMoMo(iSym,iVecVF)*(iJ-1)+(i-1)*nVir(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nVir(jSym)) = Cka(iOff1+1:iOff1+nVir(jSym))
          iOff3 = iOff3+nVir(jSym)
        end do

        do i=1,nOcc(kSym)
          iOff1 = iOffFO+nMoMo(iSym,iVecFO)*(iJ-1)+(i-1)*nFro(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nFro(jSym)) = CiK(iOff1+1:iOff1+nFro(jSym))
          iOff3 = iOff3+nFro(jSym)
          iOff1 = iOffOO+nMoMo(iSym,iVecOO)*(iJ-1)+(i-1)*nOcc(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nOcc(jSym)) = Cij(iOff1+1:iOff1+nOcc(jSym))
          iOff3 = iOff3+nOcc(jSym)
          iOff1 = iOffVO+nMoMo(iSym,iVecVO)*(iJ-1)+(i-1)*nVir(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nVir(jSym)) = Cia(iOff1+1:iOff1+nVir(jSym))
          iOff3 = iOff3+nVir(jSym)
        end do

        do i=1,nVir(kSym)
          iOff1 = iOffFV+nMoMo(iSym,iVecFV)*(iJ-1)+(i-1)*nFro(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nFro(jSym)) = CaK(iOff1+1:iOff1+nFro(jSym))
          iOff3 = iOff3+nFro(jSym)
          iOff1 = iOffOV+nMoMo(iSym,iVecOV)*(iJ-1)+(i-1)*nOcc(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nOcc(jSym)) = Cai(iOff1+1:iOff1+nOcc(jSym))
          iOff3 = iOff3+nOcc(jSym)
          iOff1 = iOffVV+nMoMo(iSym,iVecVV)*(iJ-1)+(i-1)*nVir(jSym)
          Cpq(iOff3+iOffBB+1:iOff3+iOffBB+nVir(jSym)) = Cab(iOff1+1:iOff1+nVir(jSym))
          iOff3 = iOff3+nVir(jSym)
        end do

        iOffFF = iOffFF+nFro(jSym)*nFro(kSym)
        iOffOF = iOffOF+nOcc(jSym)*nFro(kSym)
        iOffFO = iOffFO+nFro(jSym)*nOcc(kSym)
        iOffOO = iOffFO+nOcc(jSym)*nOcc(kSym)
        iOffVF = iOffVF+nVir(jSym)*nFro(kSym)
        iOffFV = iOffFV+nFro(jSym)*nVir(kSym)
        iOffVO = iOffVO+nVir(jSym)*nOcc(kSym)
        iOffOV = iOffOV+nOcc(jSym)*nVir(kSym)
        iOffVV = iOffVV+nVir(jSym)*nVir(kSym)
        iOffBB = iOffBB+nBas(jSym)*nBas(kSym)
      end do
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      do jSym=1,nSym
        kSym = ieor(iSym-1,jSym-1)+1

        ! Do backtransformations of the Cholesky-vectors to get C_mn.
        ! ----------------------------------------------------------
        !
        ! (C^-1)^T x C_pq^J = C_nq^J

        iOff1 = nLRo(iSym)*(iJ-1)+iOffLRo(iSym,jSym)
        call dGemm_('T','N',nBas(jSym),nOrb(kSym),nOrb(jSym),One,CMO_inv(1+iOffCInv(jSym)),nOrb(jSym),Cpq(1+iOff1),nOrb(jSym), &
                    Zero,Cpn,nBas(jSym))

        ! C_nq^J x  (C^-1) = C_nm^J

        iOff2 = nLRb(iSym)*(iJ-1)+iOffLRb(iSym,jSym)
        call dGemm_('N','N',nBas(jSym),nBas(kSym),nOrb(kSym),One,Cpn,nBas(jSym),CMO_inv(1+iOffCInv(kSym)),nOrb(kSym),Zero, &
                    Cmn(1+iOff2),nBas(jSym))

      end do
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      do jSym=1,nSym
        kSym = ieor(iSym-1,jSym-1)+1

        ! C_nm^J x D(MP2)_mk = U_nk^J, Eq. (43)

        iOff = nLRb(iSym)*(iJ-1)+iOffLRb(iSym,jSym)
        iOff0 = iOffD(kSym)
        call dGemm_('N','N',nBas(jSym),nBas(kSym),nBas(kSym),One,Cmn(1+iOff),nBas(jSym),MP2Density(1+iOff0),nBas(kSym),Zero, &
                    Ukn(1+iOff),nBas(jSym))

        ! D(HF)_kn x C_nm^J = V_km^J, Eq. (42)

        iOff0 = iOffD(jSym)
        call dGemm_('N','N',nBas(jSym),nBas(kSym),nBas(jSym),One,SCFDensity(1+iOff0),nBas(jSym),Cmn(1+iOff),nBas(jSym),Zero, &
                    Vkn(1+iOff),nBas(jSym))

      end do
    end do ! iJ
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Write C_nm^J, U_nk^J, and V_km^J to disk.

    lTot = nLRb(iSym)*nJ

    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuLVec,iWr,Cmn,lTot,iAdr)
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuXVec,iWr,Ukn,lTot,iAdr)
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuYVec,iWr,Vkn,lTot,iAdr)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    do lSym=1,nSym
      do iiK=1,nMP2Vec(lSym),nVec
        nK = min(nVec,nMP2Vec(lSym)-(iiK-1))

        ! Read the R-Vectors

        lTot = nMoMo(lSym,iVecOV)*nK
        iAdr = nMoMo(lSym,iVecOV)*(iiK-1)+1
        call dDaFile(lUnit_F(lSym,iTypR),iRd,Ria,lTot,iAdr)

        ! Do backtransformations of the Amplitude-vectors to get Rmnx.
        ! ------------------------------------------------------------
        ! See Eq. (39).

        iOff_Ria = 0
        iOff_Rin = 0
        iOff_Rmn = 0
        do iK=1,nK

          do jSym=1,nSym
            kSym = ieor(lSym-1,jSym-1)+1

            ! C_na x R_ai^K = R_ni^K  (Temporary storage)

            iOff1 = iOff_Ria
            iOff_Ria = iOff_Ria+nVir(jSym)*nOcc(kSym)
            iOff2 = iOff_Rin
            iOff_Rin = iOff_Rin+nBas(jSym)*nOcc(kSym)
            call dGemm_('N','N',nBas(jSym),nOcc(kSym),nVir(jSym),One,CMO_v,nBas(jSym),Ria(1+iOff1),nVir(jSym),Zero,Rin(1+iOff2), &
                        nBas(jSym))

            ! R_ni^K x C_mi^T = R_nm^K

            iOff3 = iOff_Rmn
            iOff_Rmn = iOff_Rmn+nBas(jSym)*nBas(kSym)
            call dGemm_('N','T',nBas(jSym),nBas(kSym),nOcc(kSym),One,Rin(1+iOff2),nBas(jSym),CMO_o,nBas(kSym),Zero,Rmn(1+iOff3), &
                        nBas(jSym))

          end do ! jSym
        end do   ! iK

        ! Write  R_mn^K to disk

        lTot = nLRb(iSym)*nK
        iAdr = iAdrR(iSym)+nLRb(iSym)*(iiK-1)
        call dDaFile(LuRVec,iWr,Rmn,lTot,iAdr)
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Sum_mn (R_mn^K)^T x C_mn^J = W_KJ  (Eq. 38)

        if (iSym == lSym) then
          iOffZ = iOff_WJL(iSym)+nMP2Vec(iSym)*(iiJ-1)+iiK-1
          call dGemm_('T','N',nK,nJ,nLRb(iSym),One,Rmn,nLRb(iSym),Cmn,nLRb(iSym),Zero,WJL(1+iOffZ),nMP2Vec(iSym))
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        do iJ=1,nJ
          iOff_Rin = 0
          do iK=1,nK

            do jSym=1,nSym
              kSym = ieor(iSym-1,jSym-1)+1
              mSym = ieor(lSym-1,kSym-1)+1

              iOffL = nLRb(iSym)*(iJ-1)+iOffLRb(iSym,jSym)
              iOffR = iOff_Rin
              iOff_Rin = iOff_Rin+nBas(kSym)*nOcc(mSym)

              ! C_mn^J x R_ni^K = W_mi^JK, see Eq. (44)

              call dGemm_('N','N',nBas(jSym),nOcc(mSym),nBas(kSym),One,Cmn(1+iOffL),nBas(jSym),Rin(1+iOffR),nBas(kSym),Zero,WmjKJ, &
                          nBas(jSym))

              Fac = One
              if ((iK == 1) .and. (iiK == 1)) Fac = Zero

              ! R_nm^K x W_mi^JK, last two summation in Eqs. 40 and 41

              mSym2 = ieor(lSym-1,jSym-1)+1
              iOffB = nB3(iSym)*(iJ-1)+iOffB3(iSym,mSym2)
              iOffR = nLRb(lSym)*(iK-1)+iOffLRb(lSym,mSym2)
              call dGemm_('N','N',nBas(mSym2),nOcc(mSym),nBas(jSym),One,Rmn(1+iOffR),nBas(mSym2),WmjKJ,nBas(jSym),Fac, &
                          B3jl(1+iOffB),nBas(mSym2))

            end do ! jSym

          end do ! iK
        end do   ! iJ

      end do ! iiK
    end do ! lSym

    do iJ=1,nJ

      do jSym=1,nSym
        kSym = ieor(iSym-1,jSym-1)+1

        iOffB1 = nB3(iSym)*(iJ-1)+iOffB3(iSym,jSym)
        iOffB2 = nLRb(iSym)*(iJ-1)+iOffLRb(iSym,jSym)
        iOff = iOffCMOo(kSym)

        ! Complete the 3rd RHS term in Eq. 40

        call dGemm_('N','T',nBas(jSym),nBas(kSym),nOcc(kSym),-Eight,B3jl(1+iOffB1),nBas(jSym),CMO_o(1+iOff),nBas(kSym),Zero, &
                    B3kl(1+iOffB2),nBas(jSym))

      end do ! jSym
    end do   ! iJ

    ! Write matrix to file. Note that it has the same structure as
    ! the R and L matrices.

    lTot = nLRb(iSym)*nJ
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuBTmp,iWr,B3kl,lTot,iAdr)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  end do ! iiJ
end do ! iSym

call mma_deallocate(B3jl)
call mma_deallocate(Rin)
call mma_deallocate(Cpn)
call mma_deallocate(Ria)

call mma_deallocate(Cpq)
call mma_deallocate(Cab)
call mma_deallocate(Cai)
call mma_deallocate(CaK)
call mma_deallocate(Cia)
call mma_deallocate(Cij)
call mma_deallocate(CiK)
call mma_deallocate(CKa)
call mma_deallocate(CKi)
call mma_deallocate(CJK)

call mma_deallocate(SCFDensity)

iSym = 1

lB1kl = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(B1kl,lB1kl,Label='B1kl')

lB3kl_s = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(B3kl_s,lB3kl_s,Label='B3kl_s')

lA1 = NumCho(iSym)*NumCho(iSym)
call mma_allocate(A1,lA1,Label='A1')
A1(:) = Zero

lA2 = NumCho(iSym)*NumCho(iSym)
call mma_allocate(A2,lA2,Label='A2')
A2(:) = Zero

lB2kl = nBas(iSym)*nBas(iSym)*nVec
call mma_allocate(B2kl,lB2kl,Label='B2kl')

!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Start loop over batches of J

iOff_A = 0
do iSym=1,nSym
  do iiJ=1,NumCho(iSym),nVec
    nJ = min(nVec,NumCho(iSym)-(iiJ-1))

    ! Read V_km^J

    lTot = nLRb(iSym)*nJ
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuYVec,iRd,Vkn,lTot,iAdr)

    do iiK=1,NumCho(iSym),nVec
      nK = min(nVec,NumCho(iSym)-(iiK-1))

      ! Read U_km_K

      lTot = nLRb(iSym)*nK
      iAdr = iAdrB(iSym)+nLRb(iSym)*(iiK-1)
      call dDaFile(LuXVec,iRd,Ukn,lTot,iAdr)

      ! U_km^K x V_km^J, 1st term RHS eq. 41

      iOffA = iiK-1+(iiJ-1)*NumCho(iSym)+iOff_A
      Fact = One
      call dGemm_('T','N',nK,nJ,nLRb(iSym),Fact,Ukn,nLRb(iSym),Vkn,nLRb(iSym),Zero,A1(1+iOffA),NumCho(iSym))

    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    do iiK=1,nMP2Vec(iSym),nVec
      nK = min(nVec,nMP2Vec(iSym)-(iiK-1))

      ! Read R_mn^K

      lTot = nLRb(iSym)*nK
      iAdr = iAdrR(iSym)+nLRb(iSym)*(iiK-1)
      call dDaFile(LuRVec,iRd,Rmn,lTot,iAdr)

      ! R_mn_K x W_JK, 2nd RHS term Eq. 35

      Fac = One
      if (iiK == 1) Fac = Zero
      iOffZ = iOff_WJL(iSym)+iiK-1+(iiJ-1)*nMp2Vec(iSym)
      call dGemm_('N','N',nLRb(iSym),nJ,nK,-Eight,Rmn,nLRb(iSym),WJL(1+iOffZ),nMP2Vec(iSym),Fac,B2kl,nLRb(iSym))
    end do

    ! Write to disk

    lTot = nLRb(iSym)*nJ
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuB(2),iWr,B2kl,lTot,iAdr)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read 3rd RHS term in Eq. 40

    lTot = nLRb(iSym)*nJ
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuBTmp,iRd,B3kl,lTot,iAdr)

    ! Symmetrize

    do iJ=1,nJ
      iOffL = nLRb(iSym)*(iJ-1)
      do jSym=1,nSym
        kSym = ieor(iSym-1,jSym-1)+1

        do k=1,nBas(jSym)
          do l=1,nBas(kSym)
            kl = k+nBas(jSym)*(l-1)+iOffLRb(iSym,jSym)
            kl_s = kl
            lk = l+nBas(kSym)*(k-1)+iOffLRb(iSym,kSym)
            B3kl_s(1+iOffL+kl_s-1) = (B3kl(1+iOffL+kl-1)+B3kl(1+iOffL+lk-1))*Half
          end do
        end do

      end do ! jSym
    end do

    ! Write the symmetrized 3rd term to disk

    lTot = nLRb(iSym)*nJ
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuBTmp,iWr,B3kl_s,lTot,iAdr)
    !                                                                  *
    !*******************************************************************
    !                                                                  *

    do iJ=1,nJ
      iOff = nLRb(iSym)*(iJ-1)

      ! V_kn^J x D(MP2)_kn, 2nd RHS term Eq. 40

      do jSym=1,nSym
        kSym = ieor(iSym-1,jSym-1)+1

        iOff1 = iOff+iOffLRb(iSym,jSym)
        iOff2 = iOffD(kSym)
        call dGemm_('N','N',nBas(jSym),nBas(kSym),nBas(kSym),One,Vkn(1+iOff1),nBas(jSym),MP2Density(1+iOff2),nBas(kSym),Zero, &
                    B1kl(1+iOff1),nBas(jSym))
      end do
    end do

    ! Compound 2nd and 3rd RHS term in Eq. 40.

    B1kl(1:nLRb(iSym)*nJ) = B1kl(1:nLRb(iSym)*nJ)+B3kl(1:nLRb(iSym)*nJ)

    ! Write compounded terms to disk

    lTot = nLRb(iSym)*nJ
    iAdr = iAdrB(iSym)+nLRb(iSym)*(iiJ-1)
    call dDaFile(LuB(1),iWr,B1kl,lTot,iAdr)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do ! iiJ
  iOff_A = iOff_A+NumCho(iSym)**2
end do ! iSym
!
call mma_deallocate(MP2Density)
iSym = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! Add nonseparable coulombic terms to A matrix
! --------------------------------------------
!
! -4 W_MJ x W_MK    (don't ask me why the code say 8!)
!
! Second term RHS Eq. 36, note the different order of the indices.

iOff_A2 = 0
do iSym=1,nSym
  iOff1 = iOff_WJL(iSym)
  Fact = -Eight
  call dGemm_('T','N',NumCho(iSym),NumCho(iSym),nMP2Vec(iSym),Fact,WJL(1+iOff1),nMP2Vec(iSym),WJL(1+iOff1),nMP2Vec(iSym),Zero, &
              A2(1+iOff_A2),NumCho(iSym))

  iOff_A2 = iOff_A2+NumCho(iSym)**2
end do

iSym = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! Start loop over batches of L

! Nonseparable exchange, third loop, transformation and

iOff_A = 0
do iSym=1,nSym
  do iiJ=1,NumCho(iSym),nVec
    nJ = min(nVec,NumCho(iSym)-(iiJ-1))

    lTot = nLRb(iSym)*nJ
    iAdr = iAdrB(iSym)+nlRb(iSym)*(iiJ-1)
    call dDaFile(LuBTmp,iRd,B3kl_s,lTot,iAdr)

    do iiK=1,NumCho(iSym),nVec
      nK = min(nVec,NumCho(iSym)-(iiK-1))

      lTot = nLRb(iSym)*nK
      iAdr = iAdrB(iSym)+nLRb(iSym)*(iiK-1)
      call dDaFile(LuLVec,iRd,Cmn,lTot,iAdr)

      ! B_nm^L x C_nm^K
      !
      ! This is the second term of the RHS of Eq. (41)

      iOffA = iiJ-1+(iiK-1)*NumCho(iSym)+iOff_A
      Fact = One
      call dGemm_('T','N',nJ,nK,nLRb(iSym),Fact,B3kl_s,nLRb(iSym),Cmn,nLRb(iSym),One,A1(1+iOffA),NumCho(iSym))

    end do

  end do
  iOff_A = iOff_A+NumCho(iSym)**2
end do ! iSym
!                                                                      *
!***********************************************************************
!                                                                      *
! Write the Coulomic and Exchange contributions to the 2-center
! integrals to disk. These terms could be merged to just one.

iAdr1 = 1
iAdr2 = 1
iOff = 0
do iSym=1,nSym
  lTot = NumCho(iSym)*NumCho(iSym)
  call dDaFile(LuA(1),iWr,A1(1+iOff),lTot,iAdr1)
  call dDaFile(LuA(2),iWr,A2(1+iOff),lTot,iAdr2)
  iOff = iOff+lTot
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Close some files

call DaClos(LuXVec)
call DaClos(LuYVec)
call DaClos(LuRVec)
call DaClos(LuLVec)
call DaClos(LuBTmp)
do i=1,2
  call DaClos(LuA(i))
  call DaClos(LuB(i))
end do

iClos = 2
do iSym=1,nSym
  call ChoMP2_OpenF(iClos,iTypR,iSym)
  call ChoMP2_OpenF(iClos,iTypL,iSym)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory

call mma_deallocate(Cmn)
call mma_deallocate(Rmn)
call mma_deallocate(Ukn)
call mma_deallocate(Vkn)

call mma_deallocate(B1kl)
call mma_deallocate(A1)
call mma_deallocate(B2kl)
call mma_deallocate(A2)
call mma_deallocate(B3kl)
call mma_deallocate(B3kl_s)
call mma_deallocate(WmjKJ)
call mma_deallocate(WJL)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(CMO_v)
call mma_deallocate(CMO_o)
call mma_deallocate(CMO_Inv)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ChoMP2g_GradSetup
