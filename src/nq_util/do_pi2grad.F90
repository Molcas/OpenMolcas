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
! Calling    : FZero                                                   *
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

use iSD_data
use Center_Info
use Basis_Info, only: nBas
use nq_pdft, only: lft, lGGA
use nq_Grid, only: List_G
use nq_Info

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
integer list_s(2,nlist_s), list_bas(2,nlist_s), mAO, nAOs, mGrid, nP2_ontop, nGrad_Eff, nd1mo, mRho, nCMO
real*8 D1mo(nd1mo), TabMO(mAO,mGrid,nMOs), P2_ontop(nP2_ontop,mGrid), P2_ontop_d(np2_ontop,nGrad_Eff,mGrid), CMO(nCMO)
logical ft
real*8, allocatable, dimension(:,:,:,:) :: dTabMO
real*8 RhoI(mRho,mGrid)
real*8 RhoA(mRho,mGrid)
real*8, dimension(1:mRho,1:mGrid,1:nGrad_Eff) :: dRhoI, dRhoA
integer g_eff, iGrid
real*8 TabSO(mAO,mGrid,nMOs)
real*8, dimension(mGrid*NASHT) :: P2MOCube, MOs, dMOs, MOx, MOy, MOz
integer IOff1, iOff2, iOff3, nPi, iCoordOff, iGridOff, iCoord, nBasf, nOccO, nPMO3p, iOff0, iOffF, iCoord1, iCoord2, iCoord3, &
        iCoordOff1, iCoordOff2, iCoordOff3
real*8, dimension(nPMO3p) :: P2MOCubex, P2MOCubey, P2MOCubez, dMOx, dMOy, dMOz
real*8, allocatable :: TabSO2(:)
real*8 dTabMO2(nMOs)
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

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

call FZero(P2_ontop,mGrid*nP2_ontop)
call FZero(P2_ontop_d,nP2_ontop*nGrad_Eff*mGrid)
dRhoI(1:mRho,1:mGrid,1:nGrad_Eff) = 0.0d0
dRhoA(1:mRho,1:mGrid,1:nGrad_Eff) = 0.0d0
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

call mma_Allocate(dTabMO,nP2_ontop,nMOs,nGrad_eff,mgrid,Label='dTabMO')
dTabMO(:,:,:,:) = Zero

call mma_Allocate(TabSO2,nMOs*mAO*mGrid,Label='TabSO2')

do ilist_s=1,nlist_s

  call FZero(TabSO,mAO*mGrid*nMOs)

  call mk_SOs(TabSO,mAO,mGrid,nMOs,List_s,List_Bas,nList_s,iList_s)

  call ConvertTabSO(TabSO2,TabSO,mAO,mGrid,nMOs)

  do iGrid=1,mGrid
    IGridOff = (iGrid-1)*mAO*nMOs

    do iCoord=1,3
      ICoordOff = IGridOff+(iCoord-1)*nMOs
      g_eff = list_g(iCoord,ilist_s)

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
        ICoordOff1 = IGridOff+(iCoord1-1)*nMOs
        ICoordOff2 = IGridOff+(iCoord2-1)*nMOs
        ICoordOff3 = IGridOff+(iCoord3-1)*nMOs
      else
        ICoordOff1 = 0
        ICoordOff2 = 0
        ICoordOff3 = 0
      end if

      do iIrrep=0,mIrrep-1
        nOccO = nIsh(iIrrep)+nAsh(iIrrep)
        if (nOccO == 0) cycle
        nBasF = nBas(iIrrep)

        iOffF = OffBas(iIrrep)+nFro(iIrrep)

        iOff0 = iCoordOff+OffBas(iIrrep)
        call DGEMM_('T','N',nOccO,1,nBasF,1.0d0,CMO(OffBas2(iIrrep)),nBasF,TabSO2(iOff0),nBasF,0.0d0,dTabMO2,nOccO)
        call DAXPY_(nOccO,1.0d0,dTabMO2,1,dTabMO(1,iOffF,g_eff,iGrid),nP2_ontop)

        if (lft .and. lGGA) then
          iOff1 = iCoordOff1+OffBas(iIrrep)
          iOff2 = iCoordOff2+OffBas(iIrrep)
          iOff3 = iCoordOff3+OffBas(iIrrep)

          call DGEMM_('T','N',nOccO,1,nBasF,1.0d0,CMO(OffBas2(iIrrep)),nBasF,TabSO2(iOff1),nBasF,0.0d0,dTabMO2,nOccO)
          call DAXPY_(nOccO,1.0d0,dTabMO2,1,dTabMO(2,iOffF,g_eff,iGrid),nP2_ontop)

          call DGEMM_('T','N',nOccO,1,nBasF,1.0d0,CMO(OffBas2(iIrrep)),nBasF,TabSO2(iOff2),nBasF,0.0d0,dTabMO2,nOccO)
          call DAXPY_(nOccO,1.0d0,dTabMO2,1,dTabMO(3,iOffF,g_eff,iGrid),nP2_ontop)

          call DGEMM_('T','N',nOccO,1,nBasF,1.0d0,CMO(OffBas2(iIrrep)),nBasF,TabSO2(iOff3),nBasF,0.0d0,dTabMO2,nOccO)
          call DAXPY_(nOccO,1.0d0,dTabMO2,1,dTabMO(4,iOffF,g_eff,iGrid),nP2_ontop)
        end if
      end do ! iIrrep
    end do   ! iCoord
  end do     ! iGrid
end do       ! iList_s
call mma_deAllocate(TabSO2)
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
      P2_ontop(2,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(2,iGrid)
      P2_ontop(3,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(3,iGrid)
      P2_ontop(4,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(4,iGrid)
    end if
    do g_eff=1,nGrad_eff
      P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid)+4.0d0*dRhoI(1,iGrid,g_eff)*RhoI(1,iGrid)

      if ((Functional_type == GGA_type) .and. ft) then
        !******************ADD STUFF FOR FT: HERE***************
        P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid)+4.0d0*dRhoI(2,iGrid,g_eff)*RhoI(1,iGrid)+ &
                                    8.0d0*dRhoI(1,iGrid,g_eff)*RhoI(2,iGrid)

        P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid)+4.0d0*dRhoI(3,iGrid,g_eff)*RhoI(1,iGrid)+ &
                                    8.0d0*dRhoI(1,iGrid,g_eff)*RhoI(3,iGrid)

        P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid)+4.0d0*dRhoI(4,iGrid,g_eff)*RhoI(1,iGrid)+ &
                                    8.0d0*dRhoI(1,iGrid,g_eff)*RhoI(4,iGrid)

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
!
  do iGrid=1,mGrid
    P2_ontop(1,iGrid) = P2_ontop(1,iGrid)+RhoI(1,iGrid)*RhoA(1,iGrid)
    if ((Functional_type == GGA_type) .and. ft) then
      P2_ontop(2,iGrid) = P2_ontop(2,iGrid)+2.0d0*RhoI(2,iGrid)*RhoA(1,iGrid)+2.0d0*RhoI(1,iGrid)*RhoA(2,iGrid)
      P2_ontop(3,iGrid) = P2_ontop(3,iGrid)+2.0d0*RhoI(3,iGrid)*RhoA(1,iGrid)+2.0d0*RhoI(1,iGrid)*RhoA(3,iGrid)
      P2_ontop(4,iGrid) = P2_ontop(4,iGrid)+2.0d0*RhoI(4,iGrid)*RhoA(1,iGrid)+2.0d0*RhoI(1,iGrid)*RhoA(4,iGrid)
    end if ! GGA
    do g_eff=1,nGrad_Eff
      P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid)+2.0d0*RhoI(1,iGrid)*dRhoA(1,iGrid,g_eff)+ &
                                  2.0d0*dRhoI(1,iGrid,g_eff)*RhoA(1,iGrid)
      if ((Functional_type == GGA_type) .and. ft) then

        P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid)+2.0d0*dRhoI(2,iGrid,g_eff)*RhoA(1,iGrid)+ &
                                    4.0d0*RhoI(2,iGrid)*dRhoA(1,iGrid,g_eff)+4.0d0*dRhoI(1,iGrid,g_eff)*RhoA(2,iGrid)+ &
                                    2.0d0*RhoI(1,iGrid)*dRhoA(2,iGrid,g_eff)

        P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid)+2.0d0*dRhoI(3,iGrid,g_eff)*RhoA(1,iGrid)+ &
                                    4.0d0*RhoI(3,iGrid)*dRhoA(1,iGrid,g_eff)+4.0d0*dRhoI(1,iGrid,g_eff)*RhoA(3,iGrid)+ &
                                    2.0d0*RhoI(1,iGrid)*dRhoA(3,iGrid,g_eff)

        P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid)+2.0d0*dRhoI(4,iGrid,g_eff)*RhoA(1,iGrid)+ &
                                    4.0d0*RhoI(4,iGrid)*dRhoA(1,iGrid,g_eff)+4.0d0*dRhoI(1,iGrid,g_eff)*RhoA(4,iGrid)+ &
                                    2.0d0*RhoI(1,iGrid)*dRhoA(4,iGrid,g_eff)
      end if ! GGA
    end do   ! g_eff

  end do ! loop over grid points
end if   ! if Inactive

!***********************************************************************
! Active-Active part:                                                  *
!***********************************************************************
if (NumAsh /= 0) then
  nPi = nP2_ontop

  do g_eff=1,nGrad_eff
    do iGrid=1,mGrid
      IOff1 = (iGrid-1)*NASHT
      do iIrrep=0,mIrrep-1
        IOff2 = IOff_Ash(iIrrep)+1
        IOff3 = IOff_BasAct(iIrrep)+1
        call DCopy_(nAsh(iIrrep),dTabMO(1,iOff3,g_eff,iGrid),nPi,dMOs(IOff1+IOff2),1)
        if (lft .and. lGGA) then
          call DCopy_(nAsh(iIrrep),dTabMO(2,iOff3,g_eff,iGrid),nPi,dMOx(IOff1+IOff2),1)
          call DCopy_(nAsh(iIrrep),dTabMO(3,iOff3,g_eff,iGrid),nPi,dMOy(IOff1+IOff2),1)
          call DCopy_(nAsh(iIrrep),dTabMO(4,iOff3,g_eff,iGrid),nPi,dMOz(IOff1+IOff2),1)
        end if
      end do
    end do

    do iGrid=1,mGrid
      IOff1 = (iGrid-1)*NASHT
      do IIrrep=0,mIrrep-1
        IOff2 = IOff1+iOff_Ash(IIrrep)+1
        P2_ontop_d(1,g_eff,iGrid) = P2_ontop_d(1,g_eff,iGrid)+4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCube(IOff2),1)
        if (lft .and. lGGA) then
          P2_ontop_d(2,g_eff,iGrid) = P2_ontop_d(2,g_eff,iGrid)+4.0d0*ddot_(nAsh(IIrrep),dMOx(IOff2),1,P2MOCube(IOff2),1)+ &
                                      4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCubex(IOff2),1)
          P2_ontop_d(3,g_eff,iGrid) = P2_ontop_d(3,g_eff,iGrid)+4.0d0*ddot_(nAsh(IIrrep),dMOy(IOff2),1,P2MOCube(IOff2),1)+ &
                                      4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCubey(IOff2),1)
          P2_ontop_d(4,g_eff,iGrid) = P2_ontop_d(4,g_eff,iGrid)+4.0d0*ddot_(nAsh(IIrrep),dMOz(IOff2),1,P2MOCube(IOff2),1)+ &
                                      4.0d0*ddot_(nAsh(IIrrep),dMOs(IOff2),1,P2MOCubez(IOff2),1)
        end if
      end do
    end do
  end do

  do iGrid=1,mGrid
    IOff1 = (iGrid-1)*NASHT
    do IIrrep=0,mIrrep-1
      IOff2 = IOff1+iOff_Ash(IIrrep)+1
      P2_ontop(1,iGrid) = P2_ontop(1,iGrid)+ddot_(nAsh(IIrrep),MOs(IOff2),1,P2MOCube(IOff2),1)
    end do
  end do

  if (lGGA .and. lft) then
    do iGrid=1,mGrid
      IOff1 = (iGrid-1)*NASHT
      do IIrrep=0,mIrrep-1
        IOff2 = IOff1+iOff_Ash(IIrrep)+1
        P2_ontop(2,iGrid) = P2_ontop(2,iGrid)+4.0d0*ddot_(nAsh(IIrrep),MOx(IOff2),1,P2MOCube(IOff2),1)
        P2_ontop(3,iGrid) = P2_ontop(3,iGrid)+4.0d0*ddot_(nAsh(IIrrep),MOy(IOff2),1,P2MOCube(IOff2),1)
        P2_ontop(4,iGrid) = P2_ontop(4,iGrid)+4.0d0*ddot_(nAsh(IIrrep),MOz(IOff2),1,P2MOCube(IOff2),1)
        !write(6,*) 'MOz used for dPiz'
        !write(6,*) iGrid,iOff2-iOff1
        !call RecPrt(' ',' ',MOz(iOff2),1,nAsh(iIrrep))
      end do
    end do
  end if

end if
call mma_deAllocate(dTabMO)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nAOs)

end subroutine Do_Pi2Grad
