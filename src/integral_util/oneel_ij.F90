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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine OneEl_IJ(iS,jS,iPrint,Do_PGamma,xZeta,xZI,xKappa,xPCoor,Kernel,KrnlMm,Label,lOper,nComp,CoorO,nOrdOp,iChO,iStabO, &
                    nStabO,nIC,PtChrg,nGrid,iAddPot,SOInt,l_SOInt,rFinal,nFinal,Scrtch,nScrtch,ScrSph,nScrSph,Kern,nKern)
! Thomas Bondo Pedersen and Roland Lindh, February 2011.
!
! Purpose: compute symmetry adapted one-electron integrals for
!          shell doublet iS, jS.

use Index_Functions, only: nTri_Elem1
use Real_Spherical, only: ipSph, rSph
use iSD_data, only: iSD
use Basis_Info, only: DBSC, MolWgh, Shells
use Center_Info, only: DC
use Sizes_of_Seward, only: S
use Gateway_Info, only: FNMC
use Symmetry_Info, only: ChOper, nIrrep
use rmat, only: RMat_Type_Integrals
use define_af, only: AngTp
use property_label, only: PLabel
use Molcas, only: LenIn, MxAtom
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iS, jS, iPrint, nComp, lOper(nComp), nOrdOp, iChO(nComp), iStabO(0:7), nStabO, nIC, nGrid, &
                                 iAddPot, l_SOInt, nFinal, nScrtch, nScrSph, nKern
logical(kind=iwp), intent(inout) :: Do_PGamma
real(kind=wp), intent(inout) :: xZeta(*), xZI(*)
real(kind=wp), intent(_OUT_) :: xKappa(*), xPCoor(*)
procedure(int_kernel) :: Kernel
procedure(int_mem) :: KrnlMm
character(len=8), intent(in) :: Label
real(kind=wp), intent(in) :: CoorO(3,nComp), PtChrg(nGrid)
real(kind=wp), intent(out) :: SOInt(l_SOInt), Scrtch(nScrtch), ScrSph(nScrSph)
real(kind=wp), target, intent(out) :: rFinal(nFinal), Kern(nKern)
integer(kind=iwp) :: i, iab, iAng, iAO, iAtom, iBas, iCmp, iCnt, iCnttp, iComp, iDCRR(0:7), iDCRT(0:7), ii, iiC, iIrrep, ipFnl, &
                     iPrim, ipX, ipY, ipZ, iShell, iShll, iSmLbl, iSOBlk, iStabM(0:7), iuv, jAng, jAO, jBas, jCmp, jCnt, jCnttp, &
                     jPrim, jShell, jShll, kk, l_Coord, lA0, lA1, LambdT, lB0, lB1, lDCRR, lFinal, Lmbdr, lScrSph, lScrtch, mdci, &
                     mdcj, MemAux, MemBux, MemCux, MemKer, MemKrn, mSO, nAtoms, nDCRR, nDCRT, nij, nijab, nOp(2), nOrder, NrOpr, &
                     nSO, nStabM
real(kind=wp) :: A(3), B(3), Coord(3*MxAtom), Fact, RB(3), xFactor, xMass
character(len=LenIn) :: dbas
#ifdef _GEN1INT_
logical(kind=iwp) :: DO_TRAN, NATEST
#endif
integer(kind=iwp), external :: MemSO1
logical(kind=iwp), external :: EQ

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
  if (iPrint >= 29) write(u6,*) ' nSO=',nSO
  if (nSO < 1) return
  if (l_SOInt < nSO*iBas*jBas) then
    call WarningMessage(2,'OneEl_IJ: insufficient SOInt dimension!')
    call Abend()
  end if
  SOInt(1:nSO*iBas*jBas) = Zero
  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iPrim = iSD(5,iS)
  mdci = iSD(10,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
  dbas = dc(mdci)%LblCnt(1:LenIn)
  call UpCase(dbas)
  jShll = iSD(0,jS)
  jAng = iSD(1,jS)
  jPrim = iSD(5,jS)
  mdcj = iSD(10,jS)
  jCnttp = iSD(13,jS)
  jCnt = iSD(14,jS)
  B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
  if (iPrint >= 19) then
    write(u6,*) 'interacted Ato.Fun '
    write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
  end if
  lFinal = nIC*S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  if (lFinal > nFinal) then
    call WarningMessage(2,'lFinal > nFinal')
    call Abend()
  end if
  rFinal(1:lFinal) = Zero
  call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
  call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)
  call DCR(LambdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
  if (iPrint >= 19) then
    write(u6,*)
    write(u6,*) ' g      =',nIrrep
    write(u6,*) ' u      =',dc(mdci)%nStab
    write(u6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii=0,dc(mdci)%nStab-1)
    write(u6,*) ' v      =',dc(mdcj)%nStab
    write(u6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii=0,dc(mdcj)%nStab-1)
    write(u6,*) ' LambdaR=**',LmbdR
    write(u6,*) ' r      =',nDCRR
    write(u6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii=0,nDCRR-1)
    write(u6,*) ' m      =',nStabM
    write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
  end if
  nOp(1) = NrOpr(0)
  if (nDCRR >= 1) then
    do lDCRR=0,nDCRR-1
      call OA(iDCRR(lDCRR),B,RB)
      nOp(2) = NrOpr(iDCRR(lDCRR))
      if (iPrint >= 49) write(u6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),(RB(i),i=1,3)

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
                                Shells(jShll)%Exp,Shells(iShll)%Cff_p(1,1,2),Shells(jShll)%Cff_p(1,1,2),nAtoms,Coord,nComp,rFinal, &
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
            if (btest(lOper(iComp),iIrrep)) iIC = iIC+1
          end do
        else
          !write(u6,*) "Symmetry adapt component"
          call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,rFinal,iBas,jBas,nIC,iIC,SOInt(iSOBlk),mSO,nOp)
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
  if (iPrint >= 29) write(u6,*) ' nSO=',nSO
  if (nSO < 1) return
  if (l_SOInt < nSO*iBas*jBas) then
    call WarningMessage(2,'OneEl_IJ: insufficient SOInt dimension!')
    call Abend()
  end if
  SOInt(1:nSO*iBas*jBas) = Zero

  ! Shell info

  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iPrim = iSD(5,iS)
  mdci = iSD(10,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
  dbas = dc(mdci)%LblCnt(1:LenIn)
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
    write(u6,*) 'interacted Ato.Fun '
    write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
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
    MemAux = 1+3*nTri_Elem1(la0)*nTri_Elem1(lb0+1)*nIC
    la1 = la0
    lb1 = lb0+1
    MemBux = 1+3*nTri_Elem1(la1+1)*nTri_Elem1(lb1)*nIC
    if (la1 /= 0) MemBux = MemBux+3*nTri_Elem1(la1-1)*nTri_Elem1(lb1)*nIC
    if (lb0 /= 0) then
      lb1 = lb0-1
      MemAux = MemAux+3*nTri_Elem1(la0)*nTri_Elem1(lb0-1)*nIC
      MemCux = 1+3*nTri_Elem1(la1+1)*nTri_Elem1(lb1)*nIC
      if (la1 /= 0) MemCux = MemCux+3*nTri_Elem1(la1-1)*nTri_Elem1(lb1)*nIC
    else
      MemCux = 0
    end if
    MemAux = MemAux+max(MemBux,MemCux)
    MemKer = MemKer+MemAux
  end if

  MemKrn = MemKer*iPrim*jPrim
  if (MemKrn > nKern) then
    call WarningMessage(2,'MemKrn > nKern')
    write(u6,*) 'nOrdOp,iAng,jAng=',nOrdOp,iAng,jAng
    write(u6,*) 'MemKrn=',MemKrn
    write(u6,*) 'nKern=',nKern
    call Abend()
  end if
  Kern(1:MemKrn) = Zero
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Allocate memory for the final integrals all in the
  ! primitive basis.

  lFinal = nIC*iPrim*jPrim*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  if (lFinal > nFinal) then
    call WarningMessage(2,'lFinal > nFinal')
    call Abend()
  end if
  rFinal(1:lFinal) = Zero
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scratch area for contraction step

  lScrtch = max(iPrim,jPrim)*max(iBas,jBas)*nIC*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  if (lScrtch > nScrtch) then
    call WarningMessage(2,'lScrtch > nScrtch')
    call Abend()
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scratch area for the transformation to spherical gaussians

  lScrSph = nIC*iBas*jBas*nTri_Elem1(iAng)*nTri_Elem1(jAng)
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
    write(u6,*)
    write(u6,*) ' g      =',nIrrep
    write(u6,*) ' u      =',dc(mdci)%nStab
    write(u6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii=0,dc(mdci)%nStab-1)
    write(u6,*) ' v      =',dc(mdcj)%nStab
    write(u6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii=0,dc(mdcj)%nStab-1)
    write(u6,*) ' LambdaR=**',Lmbd
    write(u6,*) ' r      =',nDCRR
    write(u6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii=0,nDCRR-1)
    write(u6,*) ' m      =',nStabM
    write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
  end if
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute normalization factor

  iuv = dc(mdci)%nStab*dc(mdcj)%nStab
  if (MolWgh == 1) then
    Fact = real(nStabO,kind=wp)/real(LambdT,kind=wp)
  else if (MolWgh == 0) then
    Fact = real(iuv*nStabO,kind=wp)/real(nIrrep**2*LambdT,kind=wp)
  else
    Fact = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nirrep*LambdT,kind=wp)
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
      if (iPrint >= 49) write(u6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),(RB(i),i=1,3)

      ! Compute kappa and P.

      call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,A,RB,xKappa,xPCoor,xZI)

      ! Compute primitive integrals. Result is ordered ij,ab.

      call Kernel(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,xZeta,xZI,xKappa,xPCoor,rFinal,iPrim*jPrim,nIC,nComp,iAng,jAng, &
                  A,RB,nOrder,Kern,MemKer,CoorO,nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
      if (iPrint >= 49) call RecPrt(' Primitive Integrals',' ',rFinal,iPrim*jPrim,nTri_Elem1(iAng)*nTri_Elem1(jAng)*nIC)

      ! Transform from primitive to contracted basis functions.

      if (iPrint >= 99) then
        call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrim,iBas)
        call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrim,jBas)
      end if

      ! Transform i,jabx to jabx,I
      kk = nTri_Elem1(iAng)*nTri_Elem1(jAng)
      call DGEMM_('T','N',jPrim*kk*nIC,iBas,iPrim,One,rFinal,iPrim,Shells(iShll)%pCff,iPrim,Zero,Scrtch,jPrim*kk*nIC)
      ! Transform j,abxI to abxI,J
      call DGEMM_('T','N',kk*nIC*iBas,jBas,jPrim,One,Scrtch,jPrim,Shells(jShll)%pCff,jPrim,Zero,ScrSph,kk*nIC*iBas)

      if (iPrint >= 99) call RecPrt(' Contracted integrals in cartesians',' ',ScrSph,kk*nIC,iBas*jBas)

      ! Transform to spherical gaussians if needed.

      if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then
        ! Result comes back as xIJAB or xIJAb
        call CarSph(ScrSph,kk,iBas*jBas*nIC,rFinal,lScrSph,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                    RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,Scrtch,iCmp*jCmp)
        call DGeTmO(Scrtch,nIC,nIC,iBas*jBas*iCmp*jCmp,rFinal,iBas*jBas*iCmp*jCmp)
      else
        ! Transpose abx,IJ back to IJ,abx
        call DGeTmO(ScrSph,kk*nIC,kk*nIC,iBas*jBas,rFinal,iBas*jBas)
      end if
      if (iPrint >= 99) call RecPrt(' Contracted integrals in Sphericals',' ',rFinal,iBas*jBas,iCmp*jCmp*nIC)

      ! Tweak here for special cases

      ipFnl = 1
      if (Label == 'P_matrix') then
        nij = iBas*jBas
        nijab = nij*iCmp*jCmp
        do iab=1,iCmp*jCmp
          ipx = ipFnl+(iab-1)*nij
          ipy = ipx+nijab
          ipz = ipy+nijab
          rFinal(ipx:ipx+nij-1) = xPCoor(1:nij)
          rFinal(ipy:ipy+nij-1) = xPCoor(nij+1:2*nij)
          rFinal(ipz:ipz+nij-1) = xPCoor(2*nij+1:3*nij)
        end do
      else if (Label == 'FMMCnX') then
        rFinal(ipFnl:ipFnl+iBas*jBas*iCmp*jCmp*nIC-1) = Half*(A(1)+RB(1))
      else if (Label == 'FMMCnY') then
        rFinal(ipFnl:ipFnl+iBas*jBas*iCmp*jCmp*nIC-1) = Half*(A(2)+RB(2))
      else if (Label == 'FMMCnZ') then
        rFinal(ipFnl:ipFnl+iBas*jBas*iCmp*jCmp*nIC-1) = Half*(A(3)+RB(3))
      else if (Label == 'Kinetic') then

        ! multiply with 1/m, where m is the mass of an electron or muon.

        xfactor = One/dbsc(iCnttp)%fMass

        ! Add the Finite Nuclear Mass Correction if activated

        if (FNMC .and. EQ(A,B) .and. (dbsc(iCnttp)%Charge /= Zero)) then
          iAtom = dbsc(iCnttp)%AtmNr
          ! Get the atom mass in au (me=1)
          xMass = dbsc(iCnttp)%CntMass
          ! Substract the electron mass to get the nuclear mass.
          xMass = xMass-real(iAtom,kind=wp)
          !write(u6,*) 'xMass=',xMass
          xfactor = xfactor+One/xMass
        end if
        !write(u6,*) 'xfactor=',xfactor
        rFinal(1:iBas*jBas*iCmp*jCmp) = xfactor*rFinal(1:iBas*jBas*iCmp*jCmp)
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
            if (btest(lOper(iComp),iIrrep)) iIC = iIC+1
          end do
        else
          call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,rFinal,iBas,jBas,nIC,iIC,SOInt(iSOBlk),mSO,nOp)
          iSOBlk = iSOBlk+mSO*iBas*jBas
        end if
      end do
    end do
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Multiply with factors due to projection operators

  if (Fact /= One) SOInt(1:nSO*iBas*jBas) = Fact*SOInt(1:nSO*iBas*jBas)
  if (iPrint >= 99) then
    write(u6,*) ' Scaling SO''s',Fact
    call RecPrt(' Accumulated SO integrals',' ',SOInt,iBas*jBas,nSO)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *

end if

return

end subroutine OneEl_IJ
