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

subroutine Do_Pi2Grad(mAO,mGrid,P2_ontop,nP2_ontop,nGrad_Eff,list_s,nlist_s,list_bas,D1mo,nd1mo,TabMO,P2_ontop_d,RhoI,RhoA,mRho, &
                      nMOs,CMO,nCMO,TabSO,ft,P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,MOs,MOx,MOy,MOz)
!***********************************************************************
!                                                                      *
! Object: Calculation P2 ontop density and its derivatives             *
!                                                                      *
! Called from: Do_batch                                                *
!                                                                      *
!    INPUT:                                                            *
!   D1mo     = one-body density matrix in MO basis                     *
!   nd1mo    = size of D1mo                                            *
!   TabMO    = MO values computed on grid                              *
!   nMOs     = number of MO basis                                      *
!   mAO      = number of derivatives of AO...                          *
!   mGrid    = number of grid points                                   *
!                                                                      *
!***********************************************************************

use Basis_Info, only: nBas
use nq_pdft, only: lft, lGGA
use nq_Grid, only: List_G
use nq_Info, only: Functional_type, GGA_type, iOff_Ash, iOff_Bas, iOff_BasAct, mBas, mIrrep, nAsh, NASHT, nFro, nIsh, OffBas, &
                   OffBas2
use Index_Functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mAO, mGrid, nP2_ontop, nGrad_Eff, nlist_s, list_s(2,nlist_s), list_bas(2,nlist_s), nd1mo, mRho, &
                                 nMOs, nCMO, nPMO3p
real(kind=wp), intent(out) :: P2_ontop(nP2_ontop,mGrid), P2_ontop_d(np2_ontop,nGrad_Eff,mGrid), TabSO(mAO,mGrid,nMOs)
real(kind=wp), intent(in) :: D1mo(nd1mo), TabMO(mAO,mGrid,nMOs), CMO(nCMO), P2MOCube(NASHT,mGrid), P2MOCubex(NASHT,nPMO3p), &
                             P2MOCubey(NASHT,nPMO3p), P2MOCubez(NASHT,nPMO3p), MOs(NASHT,mGrid), MOx(NASHT,mGrid), &
                             MOy(NASHT,mGrid), MOz(NASHT,mGrid)
real(kind=wp), intent(inout) :: RhoA(mRho,mGrid), RhoI(mRho,mGrid)
logical(kind=iwp), intent(in) :: ft
integer(kind=iwp) :: g_eff, i, i_, iCoord, iCoord1, iCoord2, iCoord3, iGrid, iIrrep, ilist_s, IOff1, iOff2, iOffF, jOffA_, jOffB_, &
                     k, k_, kIrrep, kl, l, l_, lIrrepx, nBasf, nOccO, NumAsh, NumIsh
real(kind=wp), allocatable :: dMOs(:,:), dMOx(:,:), dMOy(:,:), dMOz(:,:), dRhoA(:,:,:), dRhoI(:,:,:), dTabMO(:,:,:,:), dTabMO2(:), &
                              TabSO2(:,:,:)
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
if (nP2_ontop == 4) then
  if ((mAO /= 10) .or. (mRho /= 4)) then
    call WarningMessage(2,' Something is wrong in dim. in p2cs')
    call Abend()
  end if
else if (nP2_ontop == 6) then
  if ((mAO /= 10) .or. (mRho /= 6)) then
    call WarningMessage(2,' Something is wrong in dim. in p2cs')
    call Abend()
  end if
end if

P2_ontop(:,:) = Zero
P2_ontop_d(:,:,:) = Zero
call mma_allocate(dRhoI,mRho,mGrid,nGrad_Eff,Label='dRhoI')
call mma_allocate(dRhoA,mRho,mGrid,nGrad_Eff,Label='dRhoA')
dRhoI(:,:,:) = Zero
dRhoA(:,:,:) = Zero
jOffA_ = 0
jOffB_ = 0
do iIrrep=0,mIrrep-1
  iOff_Ash(iIrrep) = jOffA_
  iOff_Bas(iIrrep) = jOffB_
  iOff_BasAct(iIrrep) = jOffB_+nIsh(iIrrep)+nFro(iIrrep)
  jOffA_ = jOffA_+nAsh(iIrrep)
  jOffB_ = jOffB_+mBas(iIrrep)
end do
!***********************************************************************
!   P(1,...) - P_2                                                     *
!   P(2,...), P(3,...), P(4,...) - grad P_2                            *
!***********************************************************************

call mma_allocate(dTabMO,nP2_ontop,nMOs,nGrad_eff,mGrid,Label='dTabMO')
dTabMO(:,:,:,:) = Zero

call mma_allocate(TabSO2,nMOs,mAO,mGrid,Label='TabSO2')
call mma_allocate(dTabMO2,nMOs,Label='dTabMO2')

do ilist_s=1,nlist_s

  TabSO(:,:,:) = Zero

  call mk_SOs(TabSO,mAO,mGrid,nMOs,List_s,List_Bas,nList_s,iList_s)

  call ConvertTabSO(TabSO2,TabSO,mAO,mGrid,nMOs)

  do iGrid=1,mGrid

    do iCoord=1,3
      g_eff = list_g(iCoord,ilist_s)

      iCoord1 = 0
      iCoord2 = 0
      iCoord3 = 0
      if (lft .and. lGGA) then
        select case (iCoord)
          case (1)
            iCoord1 = 4
            iCoord2 = 5
            iCoord3 = 6
          case (2)
            iCoord1 = 5
            iCoord2 = 7
            iCoord3 = 8
          case (3)
            iCoord1 = 6
            iCoord2 = 8
            iCoord3 = 9
        end select
      end if

      do iIrrep=0,mIrrep-1
        nOccO = nIsh(iIrrep)+nAsh(iIrrep)
        if (nOccO == 0) cycle
        nBasF = nBas(iIrrep)

        iOffF = OffBas(iIrrep)+nFro(iIrrep)

        call DGEMM_('T','N',nOccO,1,nBasF,One,CMO(OffBas2(iIrrep):),nBasF,TabSO2(OffBas(iIrrep):,iCoord,iGrid),nBasF,Zero,dTabMO2, &
                    nOccO)
        dTabMO(1,iOffF:iOffF+nOccO-1,g_eff,iGrid) = dTabMO(1,iOffF:iOffF+nOccO-1,g_eff,iGrid)+dTabMO2(1:nOccO)

        if (lft .and. lGGA) then
          call DGEMM_('T','N',nOccO,1,nBasF,One,CMO(OffBas2(iIrrep):),nBasF,TabSO2(OffBas(iIrrep):,iCoord1,iGrid),nBasF,Zero, &
                      dTabMO2,nOccO)
          dTabMO(2,iOffF:iOffF+nOccO-1,g_eff,iGrid) = dTabMO(2,iOffF:iOffF+nOccO-1,g_eff,iGrid)+dTabMO2(1:nOccO)

          call DGEMM_('T','N',nOccO,1,nBasF,One,CMO(OffBas2(iIrrep):),nBasF,TabSO2(OffBas(iIrrep):,iCoord2,iGrid),nBasF,Zero, &
                      dTabMO2,nOccO)
          dTabMO(3,iOffF:iOffF+nOccO-1,g_eff,iGrid) = dTabMO(3,iOffF:iOffF+nOccO-1,g_eff,iGrid)+dTabMO2(1:nOccO)

          call DGEMM_('T','N',nOccO,1,nBasF,One,CMO(OffBas2(iIrrep):),nBasF,TabSO2(OffBas(iIrrep):,iCoord3,iGrid),nBasF,Zero, &
                      dTabMO2,nOccO)
          dTabMO(4,iOffF:iOffF+nOccO-1,g_eff,iGrid) = dTabMO(4,iOffF:iOffF+nOccO-1,g_eff,iGrid)+dTabMO2(1:nOccO)
        end if
      end do ! iIrrep
    end do   ! iCoord
  end do     ! iGrid
end do       ! iList_s
call mma_deallocate(TabSO2)
call mma_deallocate(dTabMO2)
!***********************************************************************
! Inactive part:                                                       *
!***********************************************************************
NumIsh = 0
NumAsh = 0
do iIrrep=0,mIrrep-1
  NumIsh = NumIsh+nISh(iIrrep)
  NumAsh = NumAsh+nAsh(iIrrep)
end do

do iGrid=1,mGrid
  do iIrrep=0,mIrrep-1
    do i_=1,nISh(iIrrep)+nFro(iIrrep)
      i = iOff_Bas(iIrrep)+i_
      RhoI(1,iGrid) = RhoI(1,iGrid)+TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
      if ((Functional_type == GGA_type) .and. ft) then
        RhoI(2,iGrid) = RhoI(2,iGrid)+TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
        RhoI(3,iGrid) = RhoI(3,iGrid)+TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
        RhoI(4,iGrid) = RhoI(4,iGrid)+TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
      end if

      ! Build dRhoI
      do g_eff=1,nGrad_eff
        dRhoI(1,iGrid,g_eff) = dRhoI(1,iGrid,g_eff)+dTabMO(1,i,g_eff,iGrid)*TabMO(1,iGrid,i) !times 2 or not?

        if ((Functional_type == GGA_type) .and. ft) then
          dRhoI(2,iGrid,g_eff) = dRhoI(2,iGrid,g_eff)+dTabMO(1,i,g_eff,iGrid)*TabMO(2,iGrid,i)+ &
                                 TabMO(1,iGrid,i)*dTabMO(2,i,g_eff,iGrid)

          dRhoI(3,iGrid,g_eff) = dRhoI(3,iGrid,g_eff)+dTabMO(1,i,g_eff,iGrid)*TabMO(3,iGrid,i)+ &
                                 TabMO(1,iGrid,i)*dTabMO(3,i,g_eff,iGrid)

          dRhoI(4,iGrid,g_eff) = dRhoI(4,iGrid,g_eff)+dTabMO(1,i,g_eff,iGrid)*TabMO(4,iGrid,i)+ &
                                 TabMO(1,iGrid,i)*dTabMO(4,i,g_eff,iGrid)

        end if ! GGA
      end do   ! g_eff

    end do     ! i_
  end do       ! iIrrep
end do         ! iGrid

if (NumIsh /= 0) then
  do iGrid=1,mGrid
    P2_ontop(1,iGrid) = RhoI(1,iGrid)*RhoI(1,iGrid)

    if ((Functional_type == GGA_type) .and. ft) then
      P2_ontop(2,iGrid) = Four*RhoI(1,iGrid)*RhoI(2,iGrid)
      P2_ontop(3,iGrid) = Four*RhoI(1,iGrid)*RhoI(3,iGrid)
      P2_ontop(4,iGrid) = Four*RhoI(1,iGrid)*RhoI(4,iGrid)
    end if
    do g_eff=1,nGrad_eff
      P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid)+Four*dRhoI(1,iGrid,g_eff)*RhoI(1,iGrid)

      if ((Functional_type == GGA_type) .and. ft) then
        !******************ADD STUFF FOR FT: HERE***************
        P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid)+Four*dRhoI(2,iGrid,g_eff)*RhoI(1,iGrid)+ &
                                    Eight*dRhoI(1,iGrid,g_eff)*RhoI(2,iGrid)

        P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid)+Four*dRhoI(3,iGrid,g_eff)*RhoI(1,iGrid)+ &
                                    Eight*dRhoI(1,iGrid,g_eff)*RhoI(3,iGrid)

        P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid)+Four*dRhoI(4,iGrid,g_eff)*RhoI(1,iGrid)+ &
                                    Eight*dRhoI(1,iGrid,g_eff)*RhoI(4,iGrid)

      end if ! GGA
    end do   ! ngrad

  end do
end if

!***********************************************************************
! Active-Inactive part:                                                *
!***********************************************************************
if ((NumIsh /= 0) .and. (NumAsh /= 0)) then
  do kIrrep=0,mIrrep-1
    do k_=1,nASh(kIrrep)
      k = k_+iOff_BasAct(kIrrep)
      do lIrrepx=0,mIrrep-1
        do l_=1,nAsh(lIrrepx)
          l = l_+iOff_BasAct(lIrrepx)
          kl = iTri(k_+iOff_Ash(kIrrep),l_+iOff_Ash(lIrrepx))
          do iGrid=1,mGrid
            RhoA(1,iGrid) = RhoA(1,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
            if ((Functional_type == GGA_type) .and. ft) then
              RhoA(2,iGrid) = RhoA(2,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
              RhoA(3,iGrid) = RhoA(3,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
              RhoA(4,iGrid) = RhoA(4,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
            end if

            do g_eff=1,nGrad_eff
              dRhoA(1,iGrid,g_eff) = dRhoA(1,iGrid,g_eff)+D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(1,iGrid,l)

              !******************ADD STUFF FOR FT: HERE***************

              if ((Functional_type == GGA_type) .and. ft) then
                dRhoA(2,iGrid,g_eff) = dRhoA(2,iGrid,g_eff)+D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(2,iGrid,l)+ &
                                       D1mo(kl)*TabMO(1,iGrid,k)*dTabMO(2,l,g_eff,iGrid)

                dRhoA(3,iGrid,g_eff) = dRhoA(3,iGrid,g_eff)+D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(3,iGrid,l)+ &
                                       D1mo(kl)*TabMO(1,iGrid,k)*dTabMO(3,l,g_eff,iGrid)

                dRhoA(4,iGrid,g_eff) = dRhoA(4,iGrid,g_eff)+D1mo(kl)*dTabMO(1,k,g_eff,iGrid)*TabMO(4,iGrid,l)+ &
                                       D1mo(kl)*TabMO(1,iGrid,k)*dTabMO(4,l,g_eff,iGrid)
              end if ! GGA

            end do
          end do     ! iGrid
        end do       ! l_
      end do         ! lIrrepx
    end do           ! k_
  end do             ! kIrrep

  do iGrid=1,mGrid
    P2_ontop(1,iGrid) = P2_ontop(1,iGrid)+RhoI(1,iGrid)*RhoA(1,iGrid)
    if ((Functional_type == GGA_type) .and. ft) then
      P2_ontop(2,iGrid) = P2_ontop(2,iGrid)+Two*RhoI(2,iGrid)*RhoA(1,iGrid)+Two*RhoI(1,iGrid)*RhoA(2,iGrid)
      P2_ontop(3,iGrid) = P2_ontop(3,iGrid)+Two*RhoI(3,iGrid)*RhoA(1,iGrid)+Two*RhoI(1,iGrid)*RhoA(3,iGrid)
      P2_ontop(4,iGrid) = P2_ontop(4,iGrid)+Two*RhoI(4,iGrid)*RhoA(1,iGrid)+Two*RhoI(1,iGrid)*RhoA(4,iGrid)
    end if ! GGA
    do g_eff=1,nGrad_Eff
      P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid)+Two*RhoI(1,iGrid)*dRhoA(1,iGrid,g_eff)+ &
                                  Two*dRhoI(1,iGrid,g_eff)*RhoA(1,iGrid)
      if ((Functional_type == GGA_type) .and. ft) then

        P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid)+Two*dRhoI(2,iGrid,g_eff)*RhoA(1,iGrid)+ &
                                    Four*RhoI(2,iGrid)*dRhoA(1,iGrid,g_eff)+Four*dRhoI(1,iGrid,g_eff)*RhoA(2,iGrid)+ &
                                    Two*RhoI(1,iGrid)*dRhoA(2,iGrid,g_eff)

        P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid)+Two*dRhoI(3,iGrid,g_eff)*RhoA(1,iGrid)+ &
                                    Four*RhoI(3,iGrid)*dRhoA(1,iGrid,g_eff)+Four*dRhoI(1,iGrid,g_eff)*RhoA(3,iGrid)+ &
                                    Two*RhoI(1,iGrid)*dRhoA(3,iGrid,g_eff)

        P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid)+Two*dRhoI(4,iGrid,g_eff)*RhoA(1,iGrid)+ &
                                    Four*RhoI(4,iGrid)*dRhoA(1,iGrid,g_eff)+Four*dRhoI(1,iGrid,g_eff)*RhoA(4,iGrid)+ &
                                    Two*RhoI(1,iGrid)*dRhoA(4,iGrid,g_eff)
      end if ! GGA
    end do   ! g_eff

  end do ! loop over grid points
end if   ! if Inactive

call mma_deallocate(dRhoI)
call mma_deallocate(dRhoA)

!***********************************************************************
! Active-Active part:                                                  *
!***********************************************************************
if (NumAsh /= 0) then

  call mma_allocate(dMOs,NASHT,mGrid,Label='dMOs')
  if (lft .and. lGGA) then
    call mma_allocate(dMOx,NASHT,mGrid,Label='dMOx')
    call mma_allocate(dMOy,NASHT,mGrid,Label='dMOy')
    call mma_allocate(dMOz,NASHT,mGrid,Label='dMOz')
  end if

  do g_eff=1,nGrad_eff
    do iGrid=1,mGrid
      do iIrrep=0,mIrrep-1
        IOff1 = IOff_Ash(iIrrep)
        IOff2 = IOff_BasAct(iIrrep)
        dMOs(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = dTabMO(1,IOff2+1:IOff2+nAsh(iIrrep),g_eff,iGrid)
        if (lft .and. lGGA) then
          dMOx(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = dTabMO(2,IOff2+1:IOff2+nAsh(iIrrep),g_eff,iGrid)
          dMOy(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = dTabMO(3,IOff2+1:IOff2+nAsh(iIrrep),g_eff,iGrid)
          dMOz(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = dTabMO(4,IOff2+1:IOff2+nAsh(iIrrep),g_eff,iGrid)
        end if
      end do
    end do

    do iGrid=1,mGrid
      do IIrrep=0,mIrrep-1
        IOff1 = iOff_Ash(IIrrep)+1
        P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid)+Four*ddot_(nAsh(IIrrep),dMOs(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)
        if (lft .and. lGGA) then
          P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid)+ &
                                      Four*ddot_(nAsh(IIrrep),dMOx(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)+ &
                                      Four*ddot_(nAsh(IIrrep),dMOs(IOff1:,iGrid),1,P2MOCubex(IOff1:,iGrid),1)
          P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid)+ &
                                      Four*ddot_(nAsh(IIrrep),dMOy(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)+ &
                                      Four*ddot_(nAsh(IIrrep),dMOs(IOff1:,iGrid),1,P2MOCubey(IOff1:,iGrid),1)
          P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid)+ &
                                      Four*ddot_(nAsh(IIrrep),dMOz(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)+ &
                                      Four*ddot_(nAsh(IIrrep),dMOs(IOff1:,iGrid),1,P2MOCubez(IOff1:,iGrid),1)
        end if
      end do
    end do
  end do

  do iGrid=1,mGrid
    do IIrrep=0,mIrrep-1
      IOff1 = iOff_Ash(IIrrep)+1
      P2_ontop(1,iGrid) = P2_ontop(1,iGrid)+ddot_(nAsh(IIrrep),MOs(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)
    end do
  end do

  if (lGGA .and. lft) then
    do iGrid=1,mGrid
      do IIrrep=0,mIrrep-1
        IOff1 = iOff_Ash(IIrrep)+1
        P2_ontop(2,iGrid) = P2_ontop(2,iGrid)+Four*ddot_(nAsh(IIrrep),MOx(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)
        P2_ontop(3,iGrid) = P2_ontop(3,iGrid)+Four*ddot_(nAsh(IIrrep),MOy(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)
        P2_ontop(4,iGrid) = P2_ontop(4,iGrid)+Four*ddot_(nAsh(IIrrep),MOz(IOff1:,iGrid),1,P2MOCube(IOff1:,iGrid),1)
        !write(u6,*) 'MOz used for dPiz'
        !write(u6,*) iGrid,iOff1
        !call RecPrt(' ',' ',MOz(iOff1:,iGrid),1,nAsh(iIrrep))
      end do
    end do
  end if

  call mma_deallocate(dMOs)
  if (lft .and. lGGA) then
    call mma_deallocate(dMOx)
    call mma_deallocate(dMOy)
    call mma_deallocate(dMOz)
  end if

end if
call mma_deallocate(dTabMO)

return

end subroutine Do_Pi2Grad
