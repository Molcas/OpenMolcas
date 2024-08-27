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
! Copyright (C) 1992,2002, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RctFld(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)
!***********************************************************************
!                                                                      *
!     Driver for RctFld_                                               *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use rctfld_module, only: lMax, MM

implicit none
integer nh1
real*8 h1(nh1), TwoHam(nh1), D(nh1), RepNuc
logical First, Dff, NonEq
integer nComp
real*8, allocatable :: Vs(:,:), QV(:,:)

nComp = (lMax+1)*(lMax+2)*(lMax+3)/6
call mma_Allocate(Vs,nComp,2,Label='Vs')
call mma_Allocate(QV,nComp,2,Label='QV')

call RctFld_Internal(MM,nComp)

call mma_deallocate(Vs)
call mma_deallocate(QV)

contains

subroutine RctFld_Internal(Q_solute,nComp)
!***********************************************************************
!                                                                      *
! Object: to apply a modification to the one-electron hamiltonian due  *
!         the reaction field. The code here is direct!                 *
!         This subroutine works only if a call to GetInf has been      *
!         prior to calling this routine.                               *
!                                                                      *
!         h1: one-electron hamiltonian to be modified. Observe that    *
!             the contribution due to the reaction field is added to   *
!             this array, i.e. it should be set prior to calling this  *
!             routine.                                                 *
!                                                                      *
!         TwoHam: dito two-electron hamiltonian.                       *
!                                                                      *
!         D:  the first order density matrix                           *
!                                                                      *
!         h1, TwoHam and D are all in the SO basis.                    *
!                                                                      *
!         Observe the energy expression for the electric field -       *
!         charge distribution interaction!                             *
!                                                                      *
!         -1/2 Sum(nl) E(tot,nl)M(tot,nl)                              *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             July '92                                                 *
!                                                                      *
!             Modified for nonequilibrum calculations January 2002 (RL)*
!***********************************************************************
# ifdef _DEBUGPRINT_
  use Basis_Info, only: nBas
  use Symmetry_Info, only: nIrrep
# endif
  use External_Centers, only: XF
  use Gateway_global, only: PrPrt
  use Gateway_Info, only: PotNuc
  use Constants, only: Half, One, Zero
  use rctfld_module, only: EPS, EPSINF, rds
  implicit none
  real*8 Origin(3)
# ifdef _DEBUGPRINT_
  character(len=72) Label
  integer lOff, iIrrep, n
# endif
  integer nComp
  real*8 Q_solute(nComp,2)

  character(len=8) Label2
  real*8 FactOp(1), E_0_NN
  integer lOper(1), ixyz, iOff, nOrdOp, iMax, ip, ix, iy, iz, iSymX, iSymY, iSymZ, iTemp, nOpr, iMltpl
  integer, external :: IrrFnc, MltLbl
  real*8, external :: DDot_
  ! Statement Function
  iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6

  lOper(1) = 1
  nOrdOp = lMax
  ! Set flag so only the diagonal blocks are computed
  Prprt = .true.
  Origin(:) = Zero

  ! Generate local multipoles in the primitive basis and accumulate to
  ! global multipoles.
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Add nuclear-nuclear contribution to the potential energy

  if (First) then

    if (NonEq) then
      call Get_dArray('RCTFLD',QV,nComp*2)
      call Get_dScalar('E_0_NN',E_0_NN)
    end if
    !1)
    ! Compute M(nuc,nl), nuclear multipole moments

    do iMax=0,lMax
      ip = 1+iOff(iMax)
      call RFNuc(Origin,Q_solute(ip,1),iMax)
    end do

    ! Add contribution from XFIELD multipoles

    ! Use Vs as temporary space, it will anyway be overwritten
    if (allocated(XF)) call XFMoment(lMax,Q_solute,Vs,nComp,Origin)

#   ifdef _DEBUGPRINT_
    call RecPrt('Nuclear Multipole Moments',' ',Q_solute(1,1),1,nComp)
#   endif

    ! Solve dielectical equation(nuclear contribution), i.e.
    ! M(nuc,nl) -> E(nuc,nl)

    call dcopy_(nComp,Q_solute(1,1),1,Vs(1,1),1)
    call AppFld(Vs(1,1),rds,Eps,lMax,EpsInf,NonEq)

#   ifdef _DEBUGPRINT_
    call RecPrt('Nuclear Electric Field',' ',Vs(1,1),1,nComp)
#   endif

    ! Vnn = Vnn - 1/2 Sum(nl) E(nuc,nl)*M(nuc,nl)

    RepNuc = PotNuc-Half*DDot_(nComp,Q_solute(1,1),1,Vs(1,1),1)

    ! Add contributions due to slow counter charges

    if (NonEq) RepNuc = RepNuc+E_0_NN
#   ifdef _DEBUGPRINT_
    write(6,*) ' RepNuc=',RepNuc
#   endif
    !2)
    ! Compute contribution to the one-electron hamiltonian

    ! hpq = hpq + Sum(nl) E(nuc,nl)*<p|M(nl)|q>

#   ifdef _DEBUGPRINT_
    write(6,*) 'h1'
    lOff = 1
    do iIrrep=0,nIrrep-1
      n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
      if (n > 0) then
        write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
        call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
        lOff = lOff+n
      end if
    end do
#   endif

    ! Add potential due to slow counter charges

    if (NonEq) then
      call dcopy_(nComp,QV(1,1),1,QV(1,2),1)
      call AppFld_NonEQ_2(QV(1,2),rds,Eps,lMax,EpsInf,NonEq)
      call DaXpY_(nComp,One,QV(1,2),1,Vs(1,1),1)
    end if

    call Drv2_RF(lOper(1),Origin,nOrdOp,Vs(1,1),lMax,h1,nh1)

#   ifdef _DEBUGPRINT_
    write(6,*) 'h1(mod)'
    lOff = 1
    do iIrrep=0,nIrrep-1
      n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
      if (n > 0) then
        write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
        call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
        lOff = lOff+n
      end if
    end do
#   endif

    ! Update h1 and RepNuc_save with respect to static contributions!

    Label2 = 'PotNuc00'
    call Put_Temp(Label2,[RepNuc],1)
    Label2 = 'h1_raw  '
    call Put_Temp(Label2,h1,nh1)

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the electronic contribution to the charge distribution.

  !3)
  ! M(el,nl) =  - Sum(p,q) Dpq <p|M(nl)|q>

  nOpr = 1
  FactOp(1) = One
  ! Reset array for storage of multipole moment expansion
  call dcopy_(nComp,[Zero],0,Q_solute(1,2),1)
  do iMltpl=1,lMax
    do ix=iMltpl,0,-1
      if (mod(ix,2) == 0) then
        iSymX = 1
      else
        ixyz = 1
        iSymX = 2**IrrFnc(ixyz)
        if (Origin(1) /= Zero) iSymX = ior(iSymX,1)
      end if
      do iy=iMltpl-ix,0,-1
        if (mod(iy,2) == 0) then
          iSymY = 1
        else
          ixyz = 2
          iSymY = 2**IrrFnc(ixyz)
          if (Origin(2) /= Zero) iSymY = ior(iSymY,1)
        end if
        iz = iMltpl-ix-iy
        if (mod(iz,2) == 0) then
          iSymZ = 1
        else
          ixyz = 4
          iSymZ = 2**IrrFnc(ixyz)
          if (Origin(3) /= Zero) iSymZ = ior(iSymZ,1)
        end if

        iTemp = MltLbl(iSymX,MltLbl(iSymY,iSymZ))
        lOper(1) = ior(lOper(1),iTemp)
      end do
    end do
  end do
# ifdef _DEBUGPRINT_
  write(6,*) '1st order density'
  lOff = 1
  do iIrrep=0,nIrrep-1
    n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
    write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
    call Triprt(Label,' ',D(lOff),nBas(iIrrep))
    lOff = lOff+n
  end do
# endif

  call Drv1_RF(FactOp,nOpr,D,nh1,Origin,lOper,Q_solute(1,2),lMax)

# ifdef _DEBUGPRINT_
  call RecPrt('Electronic Multipole Moments',' ',Q_solute(1,2),1,nComp)
# endif

  ! Solve dielectical equation(electronic contribution), i.e.
  ! M(el,nl) -> E(el,nl)

  call dcopy_(nComp,Q_solute(1,2),1,Vs(1,2),1)
  call AppFld(Vs(1,2),rds,Eps,lMax,EpsInf,NonEq)
# ifdef _DEBUGPRINT_
  call RecPrt('Electronic Electric Field',' ',Vs(1,2),1,nComp)
# endif
  !4)
  ! Compute contribution to the two-electron hamiltonian.

  ! T(D)pq = T(D)pq + Sum(nl) E(el,nl)*<p|M(nl)|q>

  call Drv2_RF(lOper(1),Origin,nOrdOp,Vs(1,2),lMax,TwoHam,nh1)

# ifdef _DEBUGPRINT_
  write(6,*) 'h1(mod)'
  lOff = 1
  do iIrrep=0,nIrrep-1
    n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
    if (n > 0) then
      write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
      call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
      lOff = lOff+n
    end if
  end do
  write(6,*) 'TwoHam(mod)'
  lOff = 1
  do iIrrep=0,nIrrep-1
    n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
    write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
    call Triprt(Label,' ',TwoHam(lOff),nBas(iIrrep))
    lOff = lOff+n
  end do
  write(6,*) ' RepNuc=',RepNuc
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Write information to be used for gradient calculations or for
  ! non-equilibrium calculations

  if (.not. NonEq) then

    ! Save total solute multipole moments and total potential
    ! of the solution.

    call dcopy_(nComp,Q_solute(1,1),1,QV(1,1),1)
    call daxpy_(nComp,One,Q_solute(1,2),1,QV(1,1),1)
    call dcopy_(nComp,Vs(1,1),1,QV(1,2),1)
    call daxpy_(nComp,One,Vs(1,2),1,QV(1,2),1)
    call Put_dArray('RCTFLD',QV,nComp*2)

    ! Compute terms to be added to RepNuc for non-equilibrium
    ! calculation.

    call dcopy_(nComp,QV(1,1),1,QV(1,2),1)
    call AppFld_NonEQ_1(QV(1,2),rds,Eps,lMax,EpsInf,NonEq)
    E_0_NN = -Half*DDot_(nComp,QV(1,1),1,QV(1,2),1)

    call dcopy_(nComp,QV(1,1),1,QV(1,2),1)
    call AppFld_NonEQ_2(QV(1,2),rds,Eps,lMax,EpsInf,NonEq)
    E_0_NN = E_0_NN+DDot_(nComp,Q_solute(1,1),1,QV(1,2),1)
    call Put_dScalar('E_0_NN',E_0_NN)

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  ! Avoid unused argument warnings
  if (.false.) call Unused_logical(Dff)

end subroutine RctFld_Internal

end subroutine RctFld
