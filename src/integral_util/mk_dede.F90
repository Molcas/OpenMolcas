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
! Copyright (C) 1990-1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine mk_DeDe(FD,nFD,mFD,ipOffD,nOffD,ipDeDe,ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,DeDe,nDeDe)
!***********************************************************************
!                                                                      *
! Object: to decontract, desymmetrize the 1st order density matrix.    *
!         The memory at this point is assumed to be large enough to do *
!         the computation in core.                                     *
!         The data is structured with respect to four indices, two (my *
!         ny or i j) refer to primitives or basis functions and two (a *
!         b) refer to the components of the cartesian or spherical     *
!         harmonic gaussians.                                          *
!                                                                      *
!         The indices are here ordered canonically!!!                  *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!             Modified for Hermite-Gauss quadrature November '90       *
!             Modified for Rys quadrature November '90                 *
!             Modified for multipole moments November '90              *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for general kernel routines January '91         *
!             Modified for nonsymmetrical operators February '91       *
!             Modified for gradients October '91                       *
!             Modified to process 1st order density matrices, Dec. '92 *
!***********************************************************************

use Index_Functions, only: iTri
use iSD_data, only: iSD
use Basis_Info, only: MolWgh, Shells
use Center_Info, only: DC
use Symmetry_Info, only: iOper, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6
#ifdef _DEBUGPRINT_
use Index_Functions, only: nTri_Elem
use define_af, only: AngTp
use Basis_Info, only: nBas
use Symmetry_Info, only: ChOper
#endif

implicit none
integer(kind=iwp), intent(in) :: nFD, mFD, nOffD, ipDeDe, MaxDe, nDeDe
real(kind=wp), intent(in) :: FD(nFD,mFD)
integer(kind=iwp), intent(out) :: ipOffD(2+mFD,nOffD), mDeDe, mIndij
logical(kind=iwp), intent(in) :: Special_NoSym, DFT_Storage
real(kind=wp), intent(inout) :: DeDe(nDeDe)
integer(kind=iwp) :: i, iAO, iAOi, iBas, iBasi, iCmp, iCmpi, iDCRR(0:7), iFD, iHigh, ij, ijCmp, ijShll, ipD00, iPrim, iPrimi, &
                     ipStart, iS, iSh, iShell, iShll, iSmLbl, iuv, j, jAO, jAOj, jBas, jBasj, jCMp, jCmpj, jOffD, jpDAO, jPrim, &
                     jPrimj, jS, jSh, jShell, jShll, lDCRR, LmbdR, mdci, mdcj, nDCRR, nOp(2), nSkal, nSO
real(kind=wp) :: FactND, Temp
real(kind=wp), allocatable :: DAO(:), DSO(:), DSOc(:), DSOp(:), Scrt(:)
integer(kind=iwp), external :: iDAMax_, MemSO1, n2Tri, NrOpr
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iIrrep, jFD
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' Differential 1st order density matrix'
iFD = 1
do iIrrep=0,nIrrep-1
  write(u6,*)
  write(u6,*) 'iIrrep=',iIrrep
  do jFD=1,mFD
    write(u6,*) 'jFD=',jFD
    write(u6,*)
    call TriPrt(' Diagonal block',' ',FD(iFD,jFD),nBas(iIrrep))
  end do
  iFD = iFD+nTri_Elem(nBas(iIrrep))
end do
#endif

mIndij = 0
jOffD = 0
ipOffD(1,:) = ipD00
ipOffD(2,:) = nIrrep
ipOffD(3,:) = MaxDe
if (mFD == 2) ipOffD(4,:) = ipD00
!                                                                      *
!***********************************************************************
!                                                                      *
call Nr_Shells(nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
! Double loop over shells. These loops decide the integral type
! i.e. (dd), (dp), (pp), etc. This is the ab index.

do iS=1,nSkal
  iShll = iSD(0,iS)
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  iAO = iSD(7,iS)
  mdci = iSD(10,iS)
  iShell = iSD(11,iS)

  do jS=1,iS
    jShll = iSD(0,jS)
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    jPrim = iSD(5,jS)
    jAO = iSD(7,jS)
    mdcj = iSD(10,jS)
    jShell = iSD(11,jS)
    ijShll = iTri(iShell,jShell)

    iSmLbl = 1
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    if (nSO == 0) cycle
    !                                                                  *
    !*******************************************************************
    !                                                                  *
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iS,jS=',iS,jS
    write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iSD(1,iS)),',',AngTp(iSD(1,jS)),') *****'
#   endif

    call mma_allocate(DAO,max(iBas*jBas,iPrim*jPrim)*iCmp*jCmp,label='DAO')

    !-------Find the DCR for A and B

    iDCRR(0:nIrrep-1) = iOper(0:nIrrep-1)
    nDCRR = nIrrep
    LmbdR = 1
#   ifdef _DEBUGPRINT_
    write(u6,'(10A)') ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
#   endif

    ! Compute normalization factor due the DCR symmetrization
    ! of the two basis functions and the operator.

    iuv = dc(mdci)%nStab*dc(mdcj)%nStab
    FactNd = real(iuv,kind=wp)/real(nIrrep*LmbdR,kind=wp)
    if (MolWgh == 1) then
      FactNd = FactNd*real(nIrrep,kind=wp)/real(iuv,kind=wp)
    else if (MolWgh == 2) then
      FactNd = sqrt(real(iuv,kind=wp))*real(nIrrep,kind=wp)/real(LmbdR,kind=wp)
    end if

    ! Allocate memory for the elements of the Fock or 1st order
    ! density matrix which are associated with the current shell
    ! pair.

    call mma_allocate(DSOp,nSO*iPrim*jPrim,label='DSOp')
    call mma_allocate(DSOc,nSO*iBas*jBas,label='DSOc')
    call mma_allocate(DSO,nSO*iPrim*jPrim,label='DSO')
    call mma_allocate(Scrt,max(iBas*jBas,iPrim*jPrim),label='Scrt')

    ! Introduce canonical order for the contracted basis

    if (iShell >= jShell) then
      iSh = iShell
      jSh = jShell
      iAOi = iAO
      jAOj = jAO
      iBasi = iBas
      jBasj = jBas
      iPrimi = iPrim
      jPrimj = jPrim
      iCmpi = iCmp
      jCmpj = jCmp
    else
      iSh = jShell
      jSh = iShell
      iAOi = jAO
      jAOj = iAO
      iBasi = jBas
      jBasj = iBas
      iPrimi = jPrim
      jPrimj = iPrim
      iCmpi = jCmp
      jCmpj = iCmp
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Loop over the density and the spin-density (optional)

    do iFD=1,mFD
#     ifdef _DEBUGPRINT_
      if (iFD == 1) then
        write(u6,*) 'Processing the density'
      else
        write(u6,*) 'Processing the spin-density'
      end if
#     endif
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Gather the elements from 1st order density / Fock matrix.

      call SOGthr(DSOc,iBasi,jBasj,nSO,FD(1,iFD),n2Tri(iSmLbl),iSmLbl,iCmpi,jCmpj,iSh,jSh,iAOi,jAOj)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Project the Fock/1st order density matrix in AO
      ! basis on to the primitive basis.

#     ifdef _DEBUGPRINT_
      call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrimi,iBasi)
      call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrimj,jBasj)
#     endif

      ! Transform IJ,AB to J,ABi
      call DGEMM_('T','T',jBasj*nSO,iPrimi,iBasi,One,DSOc,iBasi,Shells(iShll)%pCff,iPrimi,Zero,DSOp,jBasj*nSO)
      ! Transform J,ABi to AB,ij
      call DGEMM_('T','T',nSO*iPrimi,jPrimj,jBasj,One,DSOp,jBasj,Shells(jShll)%pCff,jPrimj,Zero,DSO,nSO*iPrimi)
      ! Transpose to ij,AB
      call DGeTmO(DSO,nSO,nSO,iPrimi*jPrimj,DSOp,iPrimi*jPrimj)

#     ifdef _DEBUGPRINT_
      call RecPrt(' Decontracted 1st order density/Fock matrix',' ',DSOp,iPrimi*jPrimj,nSO)
#     endif
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Loops over symmetry operations.

      nOp(1) = NrOpr(0)
      if (iFD == 1) then

        ! Store away pointer to the block of density info

        ipOffD(1,ijShll) = jOffD+ipDeDe
        ipOffD(2,ijShll) = nDCRR
        if ((nIrrep == 1) .and. Special_NoSym) then
          ipOffD(3,ijShll) = iCmp*jCmp+iPrim*jPrim+1
        else
          ipOffD(3,ijShll) = iCmp*jCmp*(iBas*jBas+1)+iPrim*jPrim+1
        end if
      else

        ! Store away pointer to the block of spin-density info

        ipOffD(4,ijShll) = jOffD+ipDeDe
      end if
      mIndij = mIndij+(nIrrep/dc(mdci)%nStab)*(nIrrep/dc(mdcj)%nStab)
#     ifdef _DEBUGPRINT_
      write(u6,*) ' ipDeDe+jOffD,nDCRR,iCmp*jCmp*iBas*jBas=',ipDeDe+jOffD,nDCRR,iCmp*jCmp*iBas*jBas
#     endif
      do lDCRR=0,nDCRR-1
        nOp(2) = NrOpr(iDCRR(lDCRR))
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Desymmetrize the 1st order density matrix(contracted).

        call Desym1(iSmLbl,iCmpi,jCmpj,iSh,jSh,iAOi,jAOj,DAO,iBasi,jBasj,DSOc,nSO,nOp,Scrt)

        ! Store away result

        Temp = Zero
        ipStart = ipDeDe+jOffD
        if (DFT_Storage) then

          ! Storage format for numerical integration

          ! D(iBas*iCmp,jBas*jCmp)

          call ResortD(DAO,DeDe(ipStart),iBas,iCmp,jBas,jCmp)
          jOffD = jOffD+iBas*iCmp*jBas*jCmp

        else

          ! Storage format for direct SCF

          ! D(iBas*jBas+1,iCmp*jCmp)

          jpDAO = 1
          do ijCmp=1,iCmp*jCmp
            if ((nIrrep /= 1) .or. (.not. Special_NoSym)) then
              DeDe(ipDeDe+jOffD:ipDeDe+jOffD+iBas*jBas-1) = DAO(jpDAO:jpDAO+iBas*jBas-1)
              jOffD = jOffD+iBas*jBas
            end if
            ! Find the largest density for this angular combination
            iHigh = iDAMax_(iBas*jBas,DAO(jpDAO),1)
            DeDe(ipDeDe+jOffD) = abs(DAO(jpDAO+iHigh-1))
            if (Temp < abs(DAO(jpDAO+iHigh-1))) Temp = abs(DAO(jpDAO+iHigh-1))
            jOffD = jOffD+1
            jpDAO = jpDAO+iBas*jBas
          end do
#         ifdef _DEBUGPRINT_
          if ((nIrrep /= 1) .or. (.not. Special_NoSym)) call RecPrt(' DAO(+AMax)',' ',DeDe(ipStart),iBas*jBas+1,iCmp*jCmp)
#         endif
        end if

        if (DFT_Storage) cycle
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Desymmetrize the 1st order density matrix(primitive).

        call Desym1(iSmLbl,iCmpi,jCmpj,iSh,jSh,iAOi,jAOj,DAO,iPrimi,jPrimj,DSOp,nSO,nOp,Scrt)

        ! Change order so it follows what is used in TwoEl

        ipStart = ipDeDe+jOffD
        do j=1,jPrimj
          do i=1,iPrimi
            ij = (j-1)*iPrimi+i
            iHigh = IDAMax_(iCmpi*jCmpj,DAO(ij),iPrimi*jPrimj)-1
            DeDe(ipDeDe+jOffD) = abs(DAO(ij+iHigh*iPrimi*jPrimj))
            jOffD = jOffD+1
          end do
        end do
        ! Find the overall largest density
        DeDe(ipDeDe+jOffD) = Temp
        jOffD = jOffD+1
#       ifdef _DEBUGPRINT_
        call RecPrt(' D,prim',' ',DeDe(ipStart),iPrimi,jPrimj)
        write(u6,*) ' Max(DAO)=',Temp
#       endif
        !                                                              *
        !***************************************************************
        !                                                              *
      end do  ! lDCRR
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do  ! iFD
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call mma_deallocate(Scrt)
    call mma_deallocate(DSO)
    call mma_deallocate(DSOc)
    call mma_deallocate(DSOp)
    call mma_deallocate(DAO)
  end do
end do
mDeDe = jOffD

if (mDeDe /= nDeDe) then
  write(u6,*) 'DeDe:  mDeDe =',mDeDe,' nDeDe =',nDeDe
  call Abend()
end if

return

end subroutine mk_DeDe
