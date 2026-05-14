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
! Copyright (C) 1990-1992,1999, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Drv2_RF(llOper,Ccoor,nOrdOp,Fldxyz,lMax,h0,nh0)
!***********************************************************************
!                                                                      *
! Object: to compute the one-electron integrals. The method employed at*
!         this point is not necessarily the fastest. However, the total*
!         time for the computation of integrals will depend on the time*
!         spent in computing the two-electron integrals.               *
!         The memory at this point is assumed to be large enough to do *
!         the computation in core.                                     *
!         The data is structured with respect to four indices, two (my *
!         ny or i j) refer to primitives or basis functions and two (a *
!         b) refer to the components of the cartesian or spherical     *
!         harmonic gaussians.                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!             Modified for Hermite-Gauss quadrature November '90       *
!             Modified for Rys quadrature November '90                 *
!             Modified for multipole moments November '90              *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for general kernel routines January  91         *
!             Modified for nonsymmetrical operators February  91       *
!             Modified to reaction field calculations July  92         *
!             Modified loop structure April  99                        *
!***********************************************************************

use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Real_Spherical, only: ipSph, rSph
use iSD_data, only: iSD
use Basis_Info, only: DBSC, MolWgh, Shells
use Center_Info, only: DC
use Gateway_global, only: PrPrt
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use define_af, only: Angtp
use Symmetry_Info, only: ChOper
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: llOper, nOrdOp, lMax, nh0
real(kind=wp), intent(in) :: Ccoor(3), Fldxyz(nTri3_Elem1(lMax))
real(kind=wp), intent(inout) :: h0(nh0)
integer(kind=iwp) :: iAng, iAO, iBas, iCmp, iCnt, iCnttp, iDCRR(0:7), iDCRRT, iDCRT(0:7), iIC, iPrim, is, iShell, iShll, iSmLbl, &
                     iStabM(0:7), iStabO(0:7), iuv, jAng, jAO, jBas, jCmp, jCnt, jCnttp, jPrim, js, jShell, jShll, kk, lDCRR, &
                     lDCRT, lFinal, LmbdR, LmbdT, mdci, mdcj, MemKer, MemKrn, mSO, nComp, nDCRR, nDCRT, nFnc, nIC, nOp(2), nOrder, &
                     NrOpr, nScr1, nScr2, nSkal, nSO, nStabM, nStabO
real(kind=wp) :: A(3), B(3), Fact
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, ii
#endif
real(kind=wp), allocatable :: Fnl(:,:), Kappa(:), Kern(:), PCoor(:,:), Scr1(:), Scr2(:), SO_Int(:), Zeta(:), ZI(:)
integer(kind=iwp), external :: MemSO1, n2Tri

#ifdef _DEBUGPRINT_
write(u6,*) ' In Drv2_RF: llOper'
write(u6,'(1X,8I5)') llOper
write(u6,*) ' In Drv2_RF: n2Tri'
write(u6,'(1X,8I5)') n2Tri(llOper)
#endif

call SOS(iStabO,nStabO,llOper)

! Auxiliary memory allocation.

call mma_allocate(Zeta,S%m2Max,Label='Zeta')
call mma_allocate(ZI,S%m2Max,Label='ZI')
call mma_allocate(Kappa,S%m2Max,Label='Kappa')
call mma_allocate(PCoor,S%m2Max,3,Label='PCoor')
!                                                                      *
!***********************************************************************
!                                                                      *
call Nr_Shells(nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
! Double loop over shells. These loops decide the integral type

do iS=1,nSkal
  iShll = iSD(0,iS)
  if (Shells(iShll)%Aux) exit
  iAng = iSD(1,iS)
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  iAO = iSD(7,iS)
  mdci = iSD(10,iS)
  iShell = iSD(11,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  do jS=1,iS
    jShll = iSD(0,jS)
    jAng = iSD(1,jS)
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    jPrim = iSD(5,jS)
    jAO = iSD(7,jS)
    mdcj = iSD(10,jS)
    jShell = iSD(11,jS)
    jCnttp = iSD(13,jS)
    jCnt = iSD(14,jS)

    iSmLbl = llOper
    if (Prprt) iSmLbl = merge(1,0,btest(iSmLbl,0))
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' nSO=',nSO
#   endif
    if (nSO == 0) cycle

#   ifdef _DEBUGPRINT_
    write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
#   endif

    ! Call kernel routine to get memory requirement. Observe, however
    ! that kernels which will use the HRR will allocate that
    ! memory internally.

    call RFMem(nOrder,MemKer,iAng,jAng,nOrdOp)
    !write(u6,*) nOrder,MemKer,iAng,jAng,nOrdOp
    MemKrn = MemKer*S%m2Max
    call mma_allocate(Kern,MemKrn,Label='Kern')

    ! Allocate memory for the final integrals all in the
    ! primitive basis.
    nComp = nTri3_Elem1(lMax)
    lFinal = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
    call mma_allocate(Fnl,lFinal,nComp+1,Label='Fnl')

    ! Scratch area for contraction step

    nScr1 = max(S%MaxPrm(iAng),S%MaxPrm(jAng))*max(S%MaxBas(iAng),S%MaxBas(jAng))*nComp*nTri_Elem1(iAng)*nTri_Elem1(jAng)
    call mma_allocate(Scr1,nScr1,Label='Scr1')

    ! Scratch area for the transformation to spherical gaussians

    nScr2 = nComp*S%MaxBas(iAng)*S%MaxBas(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
    call mma_allocate(Scr2,nScr2,Label='Scr2')

    ! At this point we can compute Zeta.
    ! This is now computed in the ij or ji order.

    call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

    ! Allocate memory for SO integrals that will be generated by
    ! this batch of AO integrals.

    call mma_allocate(SO_Int,nSO*iBas*jBas,Label='SO')
    SO_Int(:) = Zero

    ! Find the DCR for A and B

    call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)

    ! Find the stabilizer for A and B

    call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' g      =',nIrrep
    write(u6,*) ' u      =',dc(mdci)%nStab
    write(u6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii=0,dc(mdci)%nStab-1)
    write(u6,*) ' v      =',dc(mdcj)%nStab
    write(u6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii=0,dc(mdcj)%nStab-1)
    write(u6,*) ' LambdaR=',LmbdR
    write(u6,*) ' r      =',nDCRR
    write(u6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii=0,nDCRR-1)
    write(u6,*) ' m      =',nStabM
    write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    write(u6,*) ' s      =',nStabO
    write(u6,'(9A)') '(S)=',(ChOper(iStabO(ii)),ii=0,nStabO-1)
    write(u6,*) ' LambdaT=',LmbdT
    write(u6,*) ' t      =',nDCRT
    write(u6,'(9A)') '(R)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
#   endif

    ! Compute normalization factor

    iuv = dc(mdci)%nStab*dc(mdcj)%nStab
    Fact = real(iuv*nStabO,kind=wp)/real(nIrrep**2*LmbdT,kind=wp)
    if (MolWgh == 1) then
      Fact = Fact*real(nIrrep,kind=wp)**2/real(iuv,kind=wp)
    else if (MolWgh == 2) then
      Fact = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nIrrep*LmbdT,kind=wp)
    end if

    ! Loops over symmetry operations.

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),dbsc(iCnttp)%Coor(1:3,iCnt),A)
      nOp(1) = NrOpr(iDCRT(lDCRT))
#     ifdef _DEBUGPRINT_
      if (jbas < -99999) write(u6,*) 'nDCRR=',nDCRR
#     endif
      do lDCRR=0,nDCRR-1
        iDCRRT = ieor(iDCRR(lDCRR),iDCRT(lDCRT))
        call OA(iDCRRT,dbsc(jCnttp)%Coor(1:3,jCnt),B)
        nOp(2) = NrOpr(ieor(iDCRT(lDCRT),iDCRR(lDCRR)))
#       ifdef _DEBUGPRINT_
        write(u6,'(A,3(3F6.2,2X))') '***** Centers A, B, & C. *****',(A(i),i=1,3),(B(i),i=1,3),(Ccoor(i),i=1,3)
#       endif

        ! Compute kappa and P.

        call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,A,B,Kappa,PCoor,ZI)

        ! Compute primitive integrals. Result is ordered ij,ab.

        call RFInt(Zeta,Kappa,Pcoor,Fnl,iPrim*jPrim,nComp,iAng,jAng,A,B,nOrder,Kern,MemKer,Ccoor,lMax)
#       ifdef _DEBUGPRINT_
        call RecPrt(' Primitive Integrals',' ',Fnl,iPrim*jPrim*nTri_Elem1(iAng)*nTri_Elem1(jAng),nComp)
#       endif

        ! Accumulate contributions due to interaction between the
        ! electric field and the multipole moments.

        nFnc = iPrim*jPrim*nTri_Elem1(iAng)*nTri_Elem1(jAng)
        Fnl(1:nFnc,nComp+1) = Zero
        call DNaXpY(nComp,nFnc,Fldxyz,1,Fnl,1,nFnc,Fnl(1,nComp+1),1,0)
#       ifdef _DEBUGPRINT_
        call RecPrt(' Solvation integrals',' ',Fnl(1,nComp+1),iPrim*jPrim,nTri_Elem1(iAng)*nTri_Elem1(jAng))
#       endif

        ! Transform from primitive to contracted basis functions.
        ! Order of transformation is fixed. It has been shown through
        ! testing that the index order ij,ab will give a performance
        ! that is up to 20% faster than the ab,ij index order.

#       ifdef _DEBUGPRINT_
        call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrim,iBas)
        call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrim,jBas)
#       endif

        ! Transform ij,x,ab to j,xabI
        kk = nTri_Elem1(iAng)*nTri_Elem1(jAng)
        call DGEMM_('T','N',jPrim*kk,iBas,iPrim,One,Fnl(1,nComp+1),iPrim,Shells(iShll)%pCff,iPrim,Zero,Scr1,jPrim*kk)
        ! Transform j,xabI to xab,IJ
        call DGEMM_('T','N',kk*iBas,jBas,jPrim,One,Scr1,jPrim,Shells(jShll)%pCff,jPrim,Zero,Fnl(1,nComp+1),kk*iBas)

#       ifdef _DEBUGPRINT_
        call RecPrt(' Contracted integrals in cartesians',' ',Fnl(1,nComp+1),kk,iBas*jBas)
#       endif

        ! Transform to spherical gaussians if needed.

        if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

          ! Result comes back as IJxAB or IJxAb
          Scr2(1:kk*iBas*jBas) = Fnl(1:kk*iBas*jBas,nComp+1)
          call CarSph(Scr2,kk,iBas*jBas,Fnl(1,nComp+1),nScr2,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                      RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,Scr1,iCmp*jCmp)
        else
          ! Transpose back to IJ,x,ab
          call DGeTmO(Fnl(1,nComp+1),kk,kk,iBas*jBas,Scr1,iBas*jBas)
        end if
#       ifdef _DEBUGPRINT_
        call RecPrt(' Contracted Integrals in Sphericals',' ',Scr1,iBas*jBas,iCmp*jCmp)

        ! At this point accumulate the batch of integrals onto the
        ! final symmetry adapted integrals.

        call RecPrt(' Accumulated SO integrals, so far...',' ',SO_Int,iBas*jBas,nSO)
#       endif

        iSmLbl = llOper
        if (Prprt) iSmLbl = merge(1,0,btest(iSmLbl,0))
        mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
        nIC = 1
        iIC = 1
        if (mSO /= 0) &
          call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,Scr1,iBas,jBas,nIC,iIC,SO_Int,mSO,nOp)

      end do
    end do

    ! Multiply with factors due to projection operators

    if (Fact /= One) SO_Int(:) = Fact*SO_Int(:)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Scaling SO''s',Fact
    call RecPrt(' Final SO integrals',' ',SO_Int,iBas*jBas,mSO)
#   endif

    ! Accumulate contribution to the Hamiltonian.

    iSmLbl = llOper
    if (Prprt) iSmLbl = merge(1,0,btest(iSmLbl,0))
    call SOAdd(SO_Int,iBas,jBas,mSO,h0,n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)

    call mma_deallocate(SO_Int)
    call mma_deallocate(Scr2)
    call mma_deallocate(Scr1)
    call mma_deallocate(Fnl)
    call mma_deallocate(Kern)
  end do
end do

call mma_deallocate(PCoor)
call mma_deallocate(Kappa)
call mma_deallocate(ZI)
call mma_deallocate(Zeta)

return

end subroutine Drv2_RF
