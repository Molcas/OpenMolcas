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

subroutine Do_PI2(D1mo,nd1mo,TabMO,mAO,mGrid,nMOs,P2_ontop,nP2_ontop,RhoI,RhoA,mRho,P2MOCube,MOs,MOx,MOy,MOz)
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

use nq_Info, only: Functional_type, GGA_type, iOff_Ash, iOff_Bas, iOff_BasAct, mBas, mIrrep, nAsh, NASHT, nFro, nIsh
use Index_Functions, only: iTri
use Constants, only: Zero, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nd1mo, mAO, mGrid, nMOs, nP2_ontop, mRho
real(kind=wp), intent(in) :: D1mo(nd1mo), TabMO(mAO,mGrid,nMOs), P2MOCube(NASHT,mGrid), MOs(NASHT,mGrid), MOx(NASHT,mGrid), &
                             MOy(NASHT,mGrid), MOz(NASHT,mGrid)
real(kind=wp), intent(out) :: P2_ontop(nP2_ontop,mGrid)
real(kind=wp), intent(inout) :: RhoI(mRho,mGrid), RhoA(mRho,mGrid)
integer(kind=iwp) :: i, i_, iGrid, iIrrep, IOff, jOffA_, jOffB_, k, k_, kIrrep, kl, l, l_, lIrrep, NumAsh, NumIsh
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
if (nP2_ontop == 4) then
  if ((mAO /= 4) .or. (mRho /= 4)) then
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
!                                                                      *
!   P(1,...) - P_2                                                     *
!   P(2,...), P(3,...), P(4,...) - grad P_2                            *
!   Not implemented:                                                   *
!   P(5,...) - grad^2 P_2                                              *
!   P(6,...) - additional part grad^2 P_2 for CS functional            *
!   P(5) and P(6) removed                                              *
!***********************************************************************

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
    !write(u6,*) ' Symm:',iIrrep
    do i_=1,nISh(iIrrep)+nFro(iIrrep)
      i = iOff_Bas(iIrrep)+i_

      !write(u6,*) ' do_p2: Inact-Inact:',iIrrep,i,TabMO(1,iGrid,i)

      RhoI(1,iGrid) = RhoI(1,iGrid)+TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
      !write(u6,'(A15,2I3,2G15.8)') 'iGrid,i,MO,RhoI',iGrid,i,TabMO(1,iGrid,i),RhoI(1,iGrid)
      !if ((Functional_type == GGA_type) .or. Do_Grad) then
      if (Functional_type == GGA_type) then
        RhoI(2,iGrid) = RhoI(2,iGrid)+TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
        RhoI(3,iGrid) = RhoI(3,iGrid)+TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
        RhoI(4,iGrid) = RhoI(4,iGrid)+TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
      end if

    end do ! i_
  end do   ! iIrrep
end do     ! iGrid

if (NumIsh /= 0) then
  do iGrid=1,mGrid
    P2_ontop(1,iGrid) = RhoI(1,iGrid)*RhoI(1,iGrid)
    !write(u6,'(A15,I3,1G28.20)') 'iGrid,P2(1)=',iGrid,P2_ontop(1,iGrid)

    !if ((Functional_type == GGA_type) .or. Do_Grad) then
    if (Functional_type == GGA_type) then
      P2_ontop(2,iGrid) = Four*RhoI(1,iGrid)*RhoI(2,iGrid)
      !write(u6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
      P2_ontop(3,iGrid) = Four*RhoI(1,iGrid)*RhoI(3,iGrid)
      !write(u6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
      P2_ontop(4,iGrid) = Four*RhoI(1,iGrid)*RhoI(4,iGrid)
      !write(u6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
    end if
    !if ((Functional_type == LDA_type) .and. Do_Grad) then
    !  !Here I must
    !  !1. transform the 2-body density matrix to AO
    !
    !  !2. Loop over effective gradients
    !  !3. Calculate P2_ontop_d(eff_Grad,iGrid)
    !
    !end if
  end do
end if

!***********************************************************************
! Active-Inactive part:                                                *
!***********************************************************************
if ((NumIsh /= 0) .and. (NumAsh /= 0)) then
  do kIrrep=0,mIrrep-1
    do k_=1,nASh(kIrrep)
      k = k_+iOff_BasAct(kIrrep)
      do lIrrep=0,mIrrep-1
        do l_=1,nAsh(lIrrep)
          l = l_+iOff_BasAct(lIrrep)
          kl = iTri(k_+iOff_Ash(kIrrep),l_+iOff_Ash(lIrrep))
          do iGrid=1,mGrid
            RhoA(1,iGrid) = RhoA(1,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
            !write(u6,'(A35,3I3,3G15.8)') 'iGrid,k,l,D1mo(kl),Tab(k),Tab(l)=',iGrid,k,l,D1mo(kl),TabMO(1,iGrid,k),TabMO(1,iGrid,l)
            !if ((Functional_type == GGA_type) .or. Do_Grad) then
            if (Functional_type == GGA_type) then
              RhoA(2,iGrid) = RhoA(2,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
              RhoA(3,iGrid) = RhoA(3,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
              !write(u6,*) 'RhoA(4,iGrid) bf =',RhoA(4,iGrid)
              RhoA(4,iGrid) = RhoA(4,iGrid)+D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
              !write(u6,*) 'D1mo(kl),Tab(1,k),Tab(1,l)=',D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
            end if
          end do ! iGrid
        end do   ! l_
      end do     ! lIrrep
    end do       ! k_
  end do         ! kIrrep

  do iGrid=1,mGrid
    P2_ontop(1,iGrid) = P2_ontop(1,iGrid)+RhoI(1,iGrid)*RhoA(1,iGrid)
    if (Functional_type == GGA_type) then
      P2_ontop(2,iGrid) = P2_ontop(2,iGrid)+Two*RhoI(2,iGrid)*RhoA(1,iGrid)+Two*RhoI(1,iGrid)*RhoA(2,iGrid)
      !write(u6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
      P2_ontop(3,iGrid) = P2_ontop(3,iGrid)+Two*RhoI(3,iGrid)*RhoA(1,iGrid)+Two*RhoI(1,iGrid)*RhoA(3,iGrid)
      !write(u6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
      P2_ontop(4,iGrid) = P2_ontop(4,iGrid)+Two*RhoI(4,iGrid)*RhoA(1,iGrid)+Two*RhoI(1,iGrid)*RhoA(4,iGrid)
      !write(u6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
    end if
  end do ! loop over grid points
end if   ! if Inactive
!***********************************************************************
! Active-Active part:                                                  *
!***********************************************************************

if (NumAsh == 0) return

!call RecPrt('P2MOCube in do_pi2',' ',P2MOCube,NASHT,mGrid)

!call RecPrt('MOs array in do_pi2',' ',MOs,NASHT,mGrid)

do iGrid=1,mGrid
  do kIrrep=0,mIrrep-1
    IOff = iOff_Ash(kIrrep)+1
    P2_ontop(1,iGrid) = P2_ontop(1,iGrid)+ddot_(nAsh(kIrrep),MOs(IOff:,iGrid),1,P2MOCube(IOff:,iGrid),1)
  end do
end do

if (Functional_type == GGA_type) then
  do iGrid=1,mGrid
    do kIrrep=0,mIrrep-1
      IOff = iOff_Ash(kIrrep)+1
      P2_ontop(2,iGrid) = P2_ontop(2,iGrid)+Four*ddot_(nAsh(kIrrep),MOx(IOff:,iGrid),1,P2MOCube(IOff:,iGrid),1)
      P2_ontop(3,iGrid) = P2_ontop(3,iGrid)+Four*ddot_(nAsh(kIrrep),MOy(IOff:,iGrid),1,P2MOCube(IOff:,iGrid),1)
      P2_ontop(4,iGrid) = P2_ontop(4,iGrid)+Four*ddot_(nAsh(kIrrep),MOz(IOff:,iGrid),1,P2MOCube(IOff:,iGrid),1)
    end do
  end do
end if

!write(u6,*) 'On-top density new code'
!write(u6,'(10(F9.6,1X))') (P2_Ontop(1,iGrid),iGrid=1,mGrid)
!write(u6,*) (P2_Ontop(1,iGrid),iGrid=1,mGrid)

return

end subroutine Do_PI2
