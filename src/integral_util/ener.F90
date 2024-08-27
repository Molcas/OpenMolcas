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
! Copyright (C) 2000, Gunnar Karlstrom                                 *
!               2000, Roland Lindh                                     *
!***********************************************************************

subroutine ener(h1,TwoHam,D,RepNuc,nh1,First,Dff,D_Tot,Grid,nGrid_,DipMom,EField,DipEff,PolEff,Cord,MaxAto,Z_Nuc,nPolComp, &
                nAnisopol,pField,tmpF)
!***********************************************************************
!                                                                      *
!     Object: to compute the needed terms from the interaction between *
!             the charge of the QM system, the dipole moments on the   *
!             grid, and the counter charge on the boundary of the      *
!             cavity. The latter can be divided into three terms,      *
!             nuclear, electronic and Langevin dipole moment terms.    *
!                                                                      *
!     Authors: G. Karlstroem                                           *
!              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
!                                                                      *
!              and                                                     *
!                                                                      *
!              R. Lindh                                                *
!              Dept. of Chem. Phys., Univ. of Lund, Sweden.            *
!                                                                      *
!              March 2000                                              *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChBas
use Basis_Info, only: nBas
use Gateway_global, only: PrPrt
use Integral_Interfaces, only: int_kernel, int_mem, OneEl_Integrals
use Constants
use stdalloc, only: mma_allocate, mma_deallocate
use rctfld_module, only: lRFCav, TK

implicit none
integer nh1, nGrid_, MaxAto, nPolComp, nAnisopol
real*8 h1(nh1), TwoHam(nh1), D(nh1), D_tot(nh1), Grid(3,nGrid_), DipMom(3,nGrid_), EField(4,nGrid_), DipEff(nGrid_), &
       PolEff(nPolComp,nGrid_), Cord(3,MaxAto), Z_Nuc(MaxAto), pField(4,nGrid_), tmpF(4,nGrid_)
real*8 RepNuc
logical First, Dff
procedure(int_kernel) :: EFInt
procedure(int_mem) :: EFMem
logical Save_tmp, NonEq
character(len=8) Label
integer, allocatable :: ips(:), lOper(:), kOper(:)
real*8, allocatable :: C_Coor(:,:)
real*8, allocatable :: Integrals(:)
integer, external :: n2Tri
real*8 EF_Grid(3), CCoor(3)
real*8 EDip2, EInt, ESelf, ENucDip, ESimple, AGSum, dX, dY, dZ, fX, fY, fZ, PFx, PFy, PFz, FDD, FTot, FTot2, x, ag, Ex, Emx, &
       RHrmt, Sig, RepHlp, EElDip, Alfa
real*8, external :: DDOt_
integer iGrid, nOrdOp, nComp, ixyz, iSymX, iSymY, iSymZ, iSymXY, iSymXZ, iSymYZ, iSyXYZ, iSymC, iComp, ix, iy, iz, ip, iSmLbl, &
        MltLbl, nInt, iSym
integer, external :: IrrFnc
#ifdef _DEBUGPRINT_
real*8 tmp_RepNuc
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! First compute contributions due to the interaction of the QM
! charge distributions with the corresponding counter charges on
! the surface of the cavity.

! Here we have three terms
! 1) nuclear-nuclear, constant: added to RepNuc once
! 2) nuclear-electronic, constant: added to h1 once
! 3) electronic-electronic, linearly dependent of the density: TwoHam

#ifdef _DEBUGPRINT_
tmp_RepNuc = RepNuc
#endif
NonEq = .false.
if (lRFCav) call RctFld(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)
!                                                                      *
!***********************************************************************
!                                                                      *
!1) Compute interaction between the dipole moment on the grid and
!   the EF on the grid (Edip2). Eself is the self interaction term
!   of the dipole moments on the grid.
!   Both terms depend on D_tot and the modifies RepNuc each iteration.

!2) Compute interaction term between the EF of the nuclear charges
!   and the dipole moments on the grid (EnucDip). This scalar terms
!   depends also on D_tot and modifies RepNuc at each iteration.

Edip2 = Zero
Eint = Zero
Eself = Zero
Enucdip = Zero
Esimple = Zero
agsum = Zero

do iGrid=1,nGrid_
  ! Langevin dipole moment
  dx = DipMom(1,iGrid)
  dy = DipMom(2,iGrid)
  dz = DipMom(3,iGrid)
  ! Total Electric Field
  fx = EField(1,iGrid)
  fy = EField(2,iGrid)
  fz = EField(3,iGrid)
  ! Electric Field from permanent multipoles
  pfx = pField(1,iGrid)
  pfy = pField(2,iGrid)
  pfz = pField(3,iGrid)

  ! Twice the energy of Langevin dipole moments in EF

  Edip2 = Edip2-fx*dx-fy*dy-fz*dz

  fDd = fx*dx+fy*dy+fz*dz   ! f dot d

  ! Do NOT halve the contribution from permanent multipoles
  Eint = Eint-Half*(fDd+pfx*dx+pfy*dy+pfz*dz)

  Esimple = Esimple-Half*(tmpF(1,iGrid)*dx+tmpF(2,iGrid)*dy+tmpF(3,iGrid)*dz)

  ! Compute self energies of dipole moments

  ftot2 = fx*fx+fy*fy+fz*fz
  ftot = sqrt(ftot2)
  x = ftot*DipEff(iGrid)*tk
  if (x <= 0.0000001d0) then
    ag = Zero
    !alang = Zero
  else
    ex = exp(x)
    emx = One/ex
    !alang = (ex+emx)/(ex-emx)-One/x
    ag = -log((ex-emx)/(Two*x))/tk
  end if

  if (iGrid > nAnisopol) then  ! isotropic
    !uind = DipEff(iGrid)*alang+ftot*PolEff(1,iGrid)
    !Eself = Eself+ftot*uind+ag-Half*ftot**2*PolEff(1,iGrid)
    Eself = Eself+fDd-Half*ftot2*PolEff(1,iGrid)
    agsum = agsum+ag
  else ! anisotropic
    Eself = Eself+Half*(fx*dx+fy*dy+fz*dz)
  end if

  ! Compute the EF at the grid point due to the nucleus

  nOrdOp = 1
  call EFNuc(Grid(1,iGrid),Z_Nuc,Cord,MaxAto,EF_Grid,nOrdOp)

  Enucdip = Enucdip-EF_Grid(1)*dx-EF_Grid(2)*dy-EF_Grid(3)*dz

end do  ! iGrid

#ifdef _DEBUGPRINT_
write(6,*) 'Esimple             =',Esimple
write(6,*) 'EnucRctfld          =',RepNuc-tmp_RepNuc
write(6,*) 'Eself               =',Eself
write(6,*) 'Ag term             =',agsum
write(6,*) 'Electrostatic energy=',Eint
write(6,*) 'EnucDip (half)      =',Half*Enucdip
#endif

RepNuc = RepNuc+Eself+agsum+Eint+Half*Enucdip

! Note:Before PS made this code compatible with XFIE, the following
! quantity was reported as Electrostatic energy: Half*Edip2
! and: RepNuc = RepNuc + Eself + Half*Edip2 + Half*Enucdip

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute interaction term between the field of the electronic
! charge and the dipole moments. This term is modified on each
! iteration.

! Retrieve the statatic h1

Label = 'h1_raw  '
call Get_Temp(Label,h1,nh1)

! Compute the modification to h1

! Evaluate electronic EF integrals at each point and add to h1

rHrmt = One
nOrdOp = 1
nComp = 3
Sig = -One
ixyz = 1
iSymX = 2**IrrFnc(ixyz)
ixyz = 2
iSymY = 2**IrrFnc(ixyz)
ixyz = 4
iSymZ = 2**IrrFnc(ixyz)
ixyz = 3
iSymXY = 2**IrrFnc(ixyz)
ixyz = 5
iSymXZ = 2**IrrFnc(ixyz)
ixyz = 6
iSymYZ = 2**IrrFnc(ixyz)
ixyz = 7
iSyXYZ = 2**IrrFnc(ixyz)

call mma_allocate(ips,nComp,Label='ips')
call mma_allocate(lOper,nComp,Label='lOper')
call mma_allocate(kOper,nComp,Label='kOper')
call mma_allocate(C_Coor,3,nComp,Label='CCoor')

Save_tmp = PrPrt
PrPrt = .true.

RepHlp = Zero
do iGrid=1,nGrid_
  write(Label,'(A,I5)') 'EF ',iGrid
  call dcopy_(3,Grid(1,iGrid),1,CCoor,1)
  iSymC = 1
  if (Ccoor(1) /= Zero) iSymC = ior(iSymC,iSymX)
  if (Ccoor(2) /= Zero) iSymC = ior(iSymC,iSymY)
  if (Ccoor(3) /= Zero) iSymC = ior(iSymC,iSymZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero)) iSymC = ior(iSymC,iSymXY)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSymXZ)
  if ((Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSymYZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ior(iSymC,iSyXYZ)

  iComp = 0
  do ix=nOrdOp,0,-1
    do iy=nOrdOp-ix,0,-1
      iComp = iComp+1
      iz = nOrdOp-ix-iy
      ixyz = 0
      if (mod(ix,2) /= 0) ixyz = ior(ixyz,1)
      if (mod(iy,2) /= 0) ixyz = ior(ixyz,2)
      if (mod(iz,2) /= 0) ixyz = ior(ixyz,4)
      iSym = 2**IrrFnc(ixyz)
      if (Ccoor(iComp) /= Zero) iSym = ior(iSym,1)
      lOper(iComp) = MltLbl(iSymC,iSym)
      kOper(iComp) = iChBas(iComp+1)
      call dcopy_(3,Ccoor,1,C_Coor(1,iComp),1)
    end do
  end do

  call OneEl_Integrals(EFInt,EFMem,Label,ips,lOper,nComp,C_Coor,nOrdOp,rHrmt,kOper,Integrals)

  Eeldip = Zero
  do iComp=1,nComp
    ip = ips(iComp)
    iSmLbl = lOper(iComp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute properties directly from integrals

    nInt = n2Tri(iSmLbl)
    if ((nInt /= 0) .and. (abs(DipMom(iComp,iGrid)) > 1.0D-20)) then
      call CmpInt(Integrals(ip),nInt,nBas,nIrrep,iSmLbl)
      if (nInt /= nh1) then
        call WarningMessage(2,'Ener: nInt /= nh1')
        write(6,*) 'nInt=',nInt
        write(6,*) 'nh1=',nh1
        call Abend()
      end if

      ! Accumulate contribution to h1

      alfa = Sig*DipMom(iComp,iGrid)
      call DaXpY_(nInt,alfa,Integrals(ip),1,h1,1)
      EelDip = EelDip-alfa*DDot_(nh1,D_tot,1,Integrals(ip),1)

    end if

    !                                                                  *
    !*******************************************************************
    !                                                                  *

  end do ! iComp

  ! Deallocate memory for integral

  call mma_deallocate(Integrals)

  RepHlp = RepHlp+EelDip*half
  !write(6,*) 'Eeldip=',Eeldip,RepHlp

end do    ! iGrid
RepNuc = RepNuc+RepHlp
#ifdef _DEBUGPRINT_
write(6,*) 'RepHlp              =',RepHlp
#endif

PrPrt = Save_tmp

call mma_deallocate(C_Coor)
call mma_deallocate(kOper)
call mma_deallocate(lOper)
call mma_deallocate(ips)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ener
