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
! Copyright (C) 2011, Thomas Bondo Pedersen                            *
!               2011, Roland Lindh                                     *
!***********************************************************************

#include "compiler_features.h"
#ifdef _IN_MODULE_

subroutine OneEl_IJ(iS,jS,iPrint,Do_PGamma,xZeta,xZI,xKappa,xPCoor,Kernel,KrnlMm,Label,lOper,nComp,CoorO,nOrdOp,iChO,iStabO, &
                    nStabO,nIC,PtChrg,nGrid,iAddPot,SOInt,l_SOInt,final,nFinal,Scrtch,nScrtch,ScrSph,nScrSph,Kern,nKern)
! Thomas Bondo Pedersen and Roland Lindh, February 2011.
!
! Purpose: compute symmetry adapted one-electron integrals for
!          shell doublet iS, jS.

use Real_Spherical, only: ipSph, rSph
use iSD_data, only: iSD
use Basis_Info, only: DBSC, Shells, MolWgh
use Center_Info, only: DC
use Sizes_of_Seward, only: S
use Gateway_Info, only: FNMC
use Symmetry_Info, only: ChOper, nIrrep
use Constants, only: Zero, One
use rmat, only: RMat_Type_Integrals
use define_af, only: AngTp
use property_label, only: PLabel

implicit none
procedure(int_kernel) :: Kernel
procedure(int_mem) :: KrnlMm
#include "Molcas.fh"
integer iS, jS, iPrint, nComp, nOrdOp, nStabO, nIC, nGrid, iAddPot, l_SOInt
logical Do_PGamma
integer nFinal, nScrtch, nScrSph
real*8, target :: final(nFinal)
real*8 Scrtch(nScrtch), ScrSph(nScrSph), SOInt(l_SOInt)
real*8 xZeta(*), xZI(*), xKappa(*), xPCoor(*)
real*8 CoorO(3,nComp), PtChrg(nGrid)
integer lOper(nComp), iChO(nComp), iStabO(0:7)
real*8 A(3), B(3), RB(3)
character(len=8) Label
character(len=LENIN) dbas
integer nOp(2), iDCRR(0:7), iDCRT(0:7), iStabM(0:7)
integer nKern
real*8, target :: Kern(nKern)
integer, external :: MemSO1
logical, external :: EQ
real*8 Coord(3*MxAtom)
#ifdef _GEN1INT_
logical NATEST, DO_TRAN
#endif
integer i, iCmp, iBas, iAO, iShell, jCmp, jBas, jAO, jShell, nSO, iComp, iSmLbl, iShll, iAng, iPrim, mdci, iCnttp, iCnt, jShll, &
        jAng, jPrim, mdcj, jCnttp, lFinal, ii, Lmbdr, nDCRR, nStabM, lDCRR, l_Coord, nAtoms, iSOBlk, iiC, mSO, iIrrep, lA0, lA1, &
        lB0, lB1, MemAux, MemBux, MemCux, MemKer, MemKrn, lScrtch, lScrSph, iuv, LambdT, kk, ipFnl, nij, nijab, iab, ipX, ipY, &
        ipZ, jj, iAtom, nDCRT, nOrder, NrOpr, jCnt
real*8 Fact, xFactor, xMass
integer :: iTwoj(0:7) = [1,2,4,8,16,32,64,128]
! Statement function
integer ixyz, nElem
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!                                                                      *
!***********************************************************************
!                                                                      *
if (Label(1:3) == 'MAG') then
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iAO = iSD(7,iS)
  iShell = iSD(11,iS)
  jCmp = iSD(2,jS)
  jBas = iSD(3,jS)
  jAO = iSD(7,jS)
  jShell = iSD(11,jS)
  nSO = 0
  B(:) = Zero
  do iComp=1,nComp
    iSmLbl = lOper(iComp)
    nSO = nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
  end do
  if (iPrint >= 29) write(6,*) ' nSO=',nSO
  if (nSO < 1) return
  if (l_SOInt < nSO*iBas*jBas) then
    call WarningMessage(2,'OneEl_IJ: insufficient SOInt dimension!')
    call Abend()
  end if
  call dCopy_(nSO*iBas*jBas,[Zero],0,SOInt,1)
  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iPrim = iSD(5,iS)
  mdci = iSD(10,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
  dbas = dc(mdci)%LblCnt(1:LENIN)
  call UpCase(dbas)
  jShll = iSD(0,jS)
  jAng = iSD(1,jS)
  jPrim = iSD(5,jS)
  mdcj = iSD(10,jS)
  jCnttp = iSD(13,jS)
  jCnt = iSD(14,jS)
  B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
  if (iPrint >= 19) then
    write(6,*) 'interacted Ato.Fun '
    write(6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
  end if
  lFinal = nIC*S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
  if (lFinal > nFinal) then
    call WarningMessage(2,'lFinal > nFinal')
    call Abend()
  end if
  call dCopy_(lFinal,[Zero],0,final,1)
  call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
  call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)
  call DCR(LambdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
  if (iPrint >= 19) then
    write(6,*)
    write(6,*) ' g      =',nIrrep
    write(6,*) ' u      =',dc(mdci)%nStab
    write(6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii=0,dc(mdci)%nStab-1)
    write(6,*) ' v      =',dc(mdcj)%nStab
    write(6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii=0,dc(mdcj)%nStab-1)
    write(6,*) ' LambdaR=**',LmbdR
    write(6,*) ' r      =',nDCRR
    write(6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii=0,nDCRR-1)
    write(6,*) ' m      =',nStabM
    write(6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
  end if
  nOp(1) = NrOpr(0)
  if (nDCRR >= 1) then
    do lDCRR=0,nDCRR-1
      call OA(iDCRR(lDCRR),B,RB)
      nOp(2) = NrOpr(iDCRR(lDCRR))
      if (iPrint >= 49) write(6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),(RB(i),i=1,3)

      call Get_nAtoms_All(nAtoms)
      l_Coord = 3*nAtoms
      call Get_dArray('Bfn Coordinates',Coord,l_Coord)
      if (nAtoms == 2) nAtoms = nAtoms+1
#     ifdef _GEN1INT_
      NATEST = (nAtoms == 2)
      if (label(4:5) == 'PX') then
        Do_Tran = .true.
      else
        Do_Tran = .false.
      end if
      read(Label(6:),'(I3)') iatom
      call test_f90mod_sgto_mag(iShell,jShell,iCmp,jCmp,iPrim,jPrim,iAng,jAng,iPrim,jPrim,mdci,mdcj,Shells(iShll)%Exp, &
                                Shells(jShll)%Exp,Shells(iShll)%Cff_p(1,1,2),Shells(jShll)%Cff_p(1,1,2),nAtoms,Coord,nComp,final, &
                                .true.,iatom,Do_Tran)
#     else
      call WarningMessage(2,'OneEl_IJ: NO Gen1int interface available!')
      call Abend()
#     endif

      iSOBlk = 1
      iIC = 1
      do iComp=1,nComp
        iSmLbl = lOper(iComp)
        mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
        if (mSO == 0) then
          do iIrrep=0,nIrrep-1
            if (iand(lOper(iComp),iTwoj(iIrrep)) /= 0) iIC = iIC+1
          end do
        else
          !write(6,*) "Symmetry adapt component"
          call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,final,iBas,jBas,nIC,iIC,SOInt(iSOBlk),mSO,nOp)
          iSOBlk = iSOBlk+mSO*iBas*jBas
        end if
      end do
    end do
  end if

else  !  MAG Integrals

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Check memory for SO integrals that will be generated by
  ! this batch of AO integrals. Init SOInt.

  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iAO = iSD(7,iS)
  iShell = iSD(11,iS)
  jCmp = iSD(2,jS)
  jBas = iSD(3,jS)
  jAO = iSD(7,jS)
  jShell = iSD(11,jS)
  nSO = 0
  do iComp=1,nComp
    iSmLbl = lOper(iComp)
    nSO = nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
  end do
  if (iPrint >= 29) write(6,*) ' nSO=',nSO
  if (nSO < 1) return
  if (l_SOInt < nSO*iBas*jBas) then
    call WarningMessage(2,'OneEl_IJ: insufficient SOInt dimension!')
    call Abend()
  end if
  call dCopy_(nSO*iBas*jBas,[Zero],0,SOInt,1)

  ! Shell info

  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iPrim = iSD(5,iS)
  mdci = iSD(10,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
  dbas = dc(mdci)%LblCnt(1:LENIN)
  call UpCase(dbas)

  jShll = iSD(0,jS)
  jAng = iSD(1,jS)
  jPrim = iSD(5,jS)
  mdcj = iSD(10,jS)
  jCnttp = iSD(13,jS)
  jCnt = iSD(14,jS)
  B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

  ! Identify if shell doublet should be computed with special
  ! R-Matrix code.

  if ((iCnttp == jCnttp) .and. (mdcj == mdci) .and. (dbas == 'DBAS')) then
    RMat_type_integrals = .true.
    if (Do_PGamma) then
      call PGamma()
      Do_PGamma = .false.
    end if
  else
    RMat_type_integrals = .false.
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iPrint >= 19) then
    write(6,*) 'interacted Ato.Fun '
    write(6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Call kernel routine to get memory requirement.

  call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)

  ! Special additional allocation for PCC integrals

  if (PLabel /= ' ') then
    la0 = iAng
    lb0 = jAng
    MemAux = 1+3*nElem(la0)*nElem(lb0+1)*nIC
    la1 = la0
    lb1 = lb0+1
    MemBux = 1+3*nElem(la1+1)*nElem(lb1)*nIC
    if (la1 /= 0) MemBux = MemBux+3*nElem(la1-1)*nElem(lb1)*nIC
    if (lb0 /= 0) then
      lb1 = lb0-1
      MemAux = MemAux+3*nElem(la0)*nElem(lb0-1)*nIC
      MemCux = 1+3*nElem(la1+1)*nElem(lb1)*nIC
      if (la1 /= 0) MemCux = MemCux+3*nElem(la1-1)*nElem(lb1)*nIC
    else
      MemCux = 0
    end if
    MemAux = MemAux+max(MemBux,MemCux)
    MemKer = MemKer+MemAux
  end if

  MemKrn = MemKer*iPrim*jPrim
  if (MemKrn > nKern) then
    call WarningMessage(2,'MemKrn > nKern')
    write(6,*) 'nOrdOp,iAng,jAng=',nOrdOp,iAng,jAng
    write(6,*) 'MemKrn=',MemKrn
    write(6,*) 'nKern=',nKern
    call Abend()
  end if
  call dCopy_(MemKrn,[Zero],0,Kern,1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Allocate memory for the final integrals all in the
  ! primitive basis.

  lFinal = nIC*iPrim*jPrim*nElem(iAng)*nElem(jAng)
  if (lFinal > nFinal) then
    call WarningMessage(2,'lFinal > nFinal')
    call Abend()
  end if
  call dCopy_(lFinal,[Zero],0,final,1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scratch area for contraction step

  lScrtch = max(iPrim,jPrim)*max(iBas,jBas)*nIC*nElem(iAng)*nElem(jAng)
  if (lScrtch > nScrtch) then
    call WarningMessage(2,'lScrtch > nScrtch')
    call Abend()
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scratch area for the transformation to spherical gaussians

  lScrSph = nIC*iBas*jBas*nElem(iAng)*nElem(jAng)
  if (lScrSph > nScrSph) then
    call WarningMessage(2,'lScrSph > nScrSph')
    call Abend()
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! At this point we can compute Zeta.
  ! This is now computed in the ij or ji order.

  call ZXia(xZeta,xZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Find the DCR for A and B

  call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)

  ! Find the stabilizer for A and B

  call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

  call DCR(LambdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

# ifdef _DEBUGPRINT_
  if (iPrint >= 19) then
    write(6,*)
    write(6,*) ' g      =',nIrrep
    write(6,*) ' u      =',dc(mdci)%nStab
    write(6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii=0,dc(mdci)%nStab-1)
    write(6,*) ' v      =',dc(mdcj)%nStab
    write(6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii=0,dc(mdcj)%nStab-1)
    write(6,*) ' LambdaR=**',Lmbd
    write(6,*) ' r      =',nDCRR
    write(6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii=0,nDCRR-1)
    write(6,*) ' m      =',nStabM
    write(6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
  end if
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute normalization factor

  iuv = dc(mdci)%nStab*dc(mdcj)%nStab
  if (MolWgh == 1) then
    Fact = dble(nStabO)/dble(LambdT)
  else if (MolWgh == 0) then
    Fact = dble(iuv*nStabO)/dble(nIrrep**2*LambdT)
  else
    Fact = sqrt(dble(iuv))*dble(nStabO)/dble(nirrep*LambdT)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Loops over symmetry operations acting on the basis.

  nOp(1) = NrOpr(0)
  if (nDCRR >= 1) then
    do lDCRR=0,nDCRR-1
      call OA(iDCRR(lDCRR),B,RB)
      nOp(2) = NrOpr(iDCRR(lDCRR))
      if (iPrint >= 49) write(6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),(RB(i),i=1,3)

      ! Compute kappa and P.

      call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,A,RB,xKappa,xPCoor,xZI)

      ! Compute primitive integrals. Result is ordered ij,ab.

      call Kernel(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,xZeta,xZI,xKappa,xPCoor,final,iPrim*jPrim,nIC,nComp,iAng,jAng,A, &
                  RB,nOrder,Kern,MemKer,CoorO,nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
      if (iPrint >= 49) call RecPrt(' Primitive Integrals',' ',final,iPrim*jPrim,nElem(iAng)*nElem(jAng)*nIC)

      ! Transform from primitive to contracted basis functions.

      if (iPrint >= 99) then
        call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrim,iBas)
        call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrim,jBas)
      end if

      ! Transform i,jabx to jabx,I
      kk = nElem(iAng)*nElem(jAng)
      call DGEMM_('T','N',jPrim*kk*nIC,iBas,iPrim,1.0d0,final,iPrim,Shells(iShll)%pCff,iPrim,0.0d0,Scrtch,jPrim*kk*nIC)
      ! Transform j,abxI to abxI,J
      call DGEMM_('T','N',kk*nIC*iBas,jBas,jPrim,1.0d0,Scrtch,jPrim,Shells(jShll)%pCff,jPrim,0.0d0,ScrSph,kk*nIC*iBas)

      if (iPrint >= 99) call RecPrt(' Contracted integrals in cartesians',' ',ScrSph,kk*nIC,iBas*jBas)

      ! Transform to spherical gaussians if needed.

      if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then
        ! Result comes back as xIJAB or xIJAb
        call CarSph(ScrSph,kk,iBas*jBas*nIC,final,lScrSph,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                    RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,Scrtch,iCmp*jCmp)
        call DGeTmO(Scrtch,nIC,nIC,iBas*jBas*iCmp*jCmp,final,iBas*jBas*iCmp*jCmp)
      else
        ! Transpose abx,IJ back to IJ,abx
        call DGeTmO(ScrSph,kk*nIC,kk*nIC,iBas*jBas,final,iBas*jBas)
      end if
      if (iPrint >= 99) call RecPrt(' Contracted integrals in Sphericals',' ',final,iBas*jBas,iCmp*jCmp*nIC)

      ! Tweak here for special cases

      ipFnl = 1
      if (Label == 'P_matrix') then
        nij = iBas*jBas
        nijab = nij*iCmp*jCmp
        do iab=1,iCmp*jCmp
          ipx = ipFnl+(iab-1)*nij
          ipy = ipx+nijab
          ipz = ipy+nijab
          call dcopy_(nij,xPCoor(1),1,final(ipx),1)
          call dcopy_(nij,xPCoor(1+nij),1,final(ipy),1)
          call dcopy_(nij,xPCoor(1+2*nij),1,final(ipz),1)
        end do
      else if (Label == 'FMMCnX') then
        do jj=1,iBas*jBas*iCmp*jCmp*nIC
          final(ipFnl+jj-1) = (A(1)+RB(1))/2.0d0
        end do
      else if (Label == 'FMMCnY') then
        do jj=1,iBas*jBas*iCmp*jCmp*nIC
          final(ipFnl+jj-1) = (A(2)+RB(2))/2.0d0
        end do
      else if (Label == 'FMMCnZ') then
        do jj=1,iBas*jBas*iCmp*jCmp*nIC
          final(ipFnl+jj-1) = (A(3)+RB(3))/2.0d0
        end do
      else if (Label == 'Kinetic') then

        ! multiply with 1/m, where m is the mass of an electron or muon.

        xfactor = One/dbsc(iCnttp)%fMass

        ! Add the Finite Nuclear Mass Correction if activated

        if (FNMC .and. EQ(A,B) .and. (dbsc(iCnttp)%Charge /= Zero)) then
          iAtom = dbsc(iCnttp)%AtmNr
          ! Get the atom mass in au (me=1)
          xMass = dbsc(iCnttp)%CntMass
          ! Substract the electron mass to get the nuclear mass.
          xMass = xMass-dble(iAtom)
          !write(6,*) 'xMass=',xMass
          xfactor = xfactor+One/xMass
        end if
        !write(6,*) 'xfactor=',xfactor
        call DScal_(iBas*jBas*iCmp*jCmp,xfactor,final,1)
      end if

      ! At this point accumulate the batch of integrals onto the
      ! final symmetry adapted integrals.

      if (iPrint >= 99) call RecPrt(' Accumulated SO integrals, so far...',' ',SOInt,iBas*jBas,nSO)

      ! Symmetry adapt component by component

      iSOBlk = 1
      iIC = 1
      do iComp=1,nComp
        iSmLbl = lOper(iComp)
        mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
        if (mSO == 0) then
          do iIrrep=0,nIrrep-1
            if (iand(lOper(iComp),iTwoj(iIrrep)) /= 0) iIC = iIC+1
          end do
        else
          call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,final,iBas,jBas,nIC,iIC,SOInt(iSOBlk),mSO,nOp)
          iSOBlk = iSOBlk+mSO*iBas*jBas
        end if
      end do
    end do
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Multiply with factors due to projection operators

  if (Fact /= One) call dScal_(nSO*iBas*jBas,Fact,SOInt,1)
  if (iPrint >= 99) then
    write(6,*) ' Scaling SO''s',Fact
    call RecPrt(' Accumulated SO integrals',' ',SOInt,iBas*jBas,nSO)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *

end if

return

end subroutine OneEl_IJ

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(OneEl_IJ)

#endif
