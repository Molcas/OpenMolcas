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

use Index_Functions, only: nTri3_Elem, nTri3_Elem1
use External_Centers, only: XF
use Gateway_global, only: PrPrt
use Gateway_Info, only: PotNuc
use rctfld_module, only: EPS, EPSINF, lMax, MM, rds
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Zero, Half
use Definitions, only: wp, iwp
# ifdef _DEBUGPRINT_
use Index_Functions, only: nTri_Elem
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Definitions, only: u6
# endif

implicit none
integer(kind=iwp), intent(in) :: nh1
real(kind=wp), intent(inout) :: h1(nh1), TwoHam(nh1), RepNuc
real(kind=wp), intent(in) :: D(nh1)
logical(kind=iwp), intent(in) :: First, Dff, NonEq
integer(kind=iwp) :: iMax, iMltpl, ip, iSymX, iSymY, iSymZ, iTemp, ix, ixyz, iy, iz, lOper(1), nComp, nOpr, nOrdOp
real(kind=wp) :: E_0_NN, FactOp(1), Origin(3)
character(len=8) :: Label2
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iIrrep, lOff, n
character(len=72) :: Label
#endif
real(kind=wp), allocatable :: QV(:,:), Vs(:,:)
integer(kind=iwp), external :: IrrFnc, MltLbl
real(kind=wp), external :: DDot_

#include "macros.fh"
unused_var(Dff)

nComp = nTri3_Elem1(lMax)
call mma_Allocate(Vs,nComp,2,Label='Vs')
call mma_Allocate(QV,nComp,2,Label='QV')

lOper(1) = 1
nOrdOp = lMax
! Set flag so only the diagonal blocks are computed
Prprt = .true.
Origin(:) = Zero

! Generate local multipoles in the primitive basis and accumulate to
! global multipoles.
!                                                                      *
!***********************************************************************
!                                                                      *
! Add nuclear-nuclear contribution to the potential energy

if (First) then

  if (NonEq) then
    call Get_dArray('RCTFLD',QV,nComp*2)
    call Get_dScalar('E_0_NN',E_0_NN)
  end if
  !1)
  ! Compute M(nuc,nl), nuclear multipole moments

  do iMax=0,lMax
    ip = 1+nTri3_Elem(iMax)
    call RFNuc(Origin,MM(ip,1),iMax)
  end do

  ! Add contribution from XFIELD multipoles

  ! Use Vs as temporary space, it will anyway be overwritten
  if (allocated(XF)) call XFMoment(lMax,MM,Vs,nComp,Origin)

# ifdef _DEBUGPRINT_
  call RecPrt('Nuclear Multipole Moments',' ',MM(:,1),1,nComp)
# endif

  ! Solve dielectical equation(nuclear contribution), i.e.
  ! M(nuc,nl) -> E(nuc,nl)

  Vs(:,1) = MM(:,1)
  call AppFld(Vs(:,1),rds,Eps,lMax,EpsInf,NonEq)

# ifdef _DEBUGPRINT_
  call RecPrt('Nuclear Electric Field',' ',Vs(:,1),1,nComp)
# endif

  ! Vnn = Vnn - 1/2 Sum(nl) E(nuc,nl)*M(nuc,nl)

  RepNuc = PotNuc-Half*DDot_(nComp,MM(:,1),1,Vs(:,1),1)

  ! Add contributions due to slow counter charges

  if (NonEq) RepNuc = RepNuc+E_0_NN
# ifdef _DEBUGPRINT_
  write(u6,*) ' RepNuc=',RepNuc
# endif
  !2)
  ! Compute contribution to the one-electron hamiltonian

  ! hpq = hpq + Sum(nl) E(nuc,nl)*<p|M(nl)|q>

# ifdef _DEBUGPRINT_
  write(u6,*) 'h1'
  lOff = 1
  do iIrrep=0,nIrrep-1
    n = nTri_Elem(nBas(iIrrep))
    if (n > 0) then
      write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
      call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
      lOff = lOff+n
    end if
  end do
# endif

  ! Add potential due to slow counter charges

  if (NonEq) then
    QV(:,2) = QV(:,1)
    call AppFld_NonEQ_2(QV(:,2),rds,Eps,lMax,EpsInf)
    Vs(:,1) = Vs(:,1)+QV(:,2)
  end if

  call Drv2_RF(lOper(1),Origin,nOrdOp,Vs(:,1),lMax,h1,nh1)

# ifdef _DEBUGPRINT_
  write(u6,*) 'h1(mod)'
  lOff = 1
  do iIrrep=0,nIrrep-1
    n = nTri_Elem(nBas(iIrrep))
    if (n > 0) then
      write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
      call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
      lOff = lOff+n
    end if
  end do
# endif

  ! Update h1 and RepNuc_save with respect to static contributions!

  Label2 = 'PotNuc00'
  call Put_Temp(Label2,[RepNuc],1)
  Label2 = 'h1_raw  '
  call Put_Temp(Label2,h1,nh1)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the electronic contribution to the charge distribution.

!3)
! M(el,nl) =  - Sum(p,q) Dpq <p|M(nl)|q>

nOpr = 1
FactOp(1) = One
! Reset array for storage of multipole moment expansion
MM(:,2) = Zero
do iMltpl=1,lMax
  do ix=iMltpl,0,-1
  if (mod(ix,2) == 0) then
        iSymX = 1
    else
      ixyz = 1
      iSymX = 2**IrrFnc(ixyz)
      if (Origin(1) /= Zero) iSymX = ibset(iSymX,0)
    end if
    do iy=iMltpl-ix,0,-1
      if (mod(iy,2) == 0) then
        iSymY = 1
      else
        ixyz = 2
        iSymY = 2**IrrFnc(ixyz)
        if (Origin(2) /= Zero) iSymY = ibset(iSymY,0)
      end if
      iz = iMltpl-ix-iy
      if (mod(iz,2) == 0) then
        iSymZ = 1
      else
        ixyz = 4
        iSymZ = 2**IrrFnc(ixyz)
        if (Origin(3) /= Zero) iSymZ = ibset(iSymZ,0)
      end if

      iTemp = MltLbl(iSymX,MltLbl(iSymY,iSymZ))
      lOper(1) = ior(lOper(1),iTemp)
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) '1st order density'
lOff = 1
do iIrrep=0,nIrrep-1
  n = nTri_Elem(nBas(iIrrep))
  write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
  call Triprt(Label,' ',D(lOff),nBas(iIrrep))
  lOff = lOff+n
end do
#endif

call Drv1_RF(FactOp,nOpr,D,nh1,Origin,lOper,MM(:,2),lMax)

#ifdef _DEBUGPRINT_
call RecPrt('Electronic Multipole Moments',' ',MM(:,2),1,nComp)
#endif

! Solve dielectical equation(electronic contribution), i.e.
! M(el,nl) -> E(el,nl)

Vs(:,2) = MM(:,2)
call AppFld(Vs(:,2),rds,Eps,lMax,EpsInf,NonEq)
#ifdef _DEBUGPRINT_
call RecPrt('Electronic Electric Field',' ',Vs(:,2),1,nComp)
#endif
!4)
! Compute contribution to the two-electron hamiltonian.

! T(D)pq = T(D)pq + Sum(nl) E(el,nl)*<p|M(nl)|q>

call Drv2_RF(lOper(1),Origin,nOrdOp,Vs(:,2),lMax,TwoHam,nh1)

#ifdef _DEBUGPRINT_
write(u6,*) 'h1(mod)'
lOff = 1
do iIrrep=0,nIrrep-1
  n = nTri_Elem(nBas(iIrrep))
  if (n > 0) then
    write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
    call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
    lOff = lOff+n
  end if
end do
write(u6,*) 'TwoHam(mod)'
lOff = 1
do iIrrep=0,nIrrep-1
  n = nTri_Elem(nBas(iIrrep))
  write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
  call Triprt(Label,' ',TwoHam(lOff),nBas(iIrrep))
  lOff = lOff+n
end do
write(u6,*) ' RepNuc=',RepNuc
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Write information to be used for gradient calculations or for
! non-equilibrium calculations

if (.not. NonEq) then

  ! Save total solute multipole moments and total potential
  ! of the solution.

  QV(:,1) = MM(:,1)+MM(:,2)
  QV(:,2) = Vs(:,1)+Vs(:,2)
  call Put_dArray('RCTFLD',QV,nComp*2)

  ! Compute terms to be added to RepNuc for non-equilibrium
  ! calculation.

  QV(:,2) = QV(:,1)
  call AppFld_NonEQ_1(QV(:,2),rds,Eps,lMax,EpsInf)
  E_0_NN = -Half*DDot_(nComp,QV(:,1),1,QV(:,2),1)

  QV(:,2) = QV(:,1)
  call AppFld_NonEQ_2(QV(:,2),rds,Eps,lMax,EpsInf)
  E_0_NN = E_0_NN+DDot_(nComp,MM(:,1),1,QV(:,2),1)
  call Put_dScalar('E_0_NN',E_0_NN)

end if
!                                                                      *
!***********************************************************************
!                                                                      *

call mma_deallocate(Vs)
call mma_deallocate(QV)

end subroutine RctFld
