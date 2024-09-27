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

!#define _DEBUGPRINT_
subroutine eperm(D_Tot,nDens,Ravxyz,Cavxyz,nCavxyz_,dEF,Grid,nGrid_,Cord,MaxAto,Z_Nuc,xfEF)
!***********************************************************************
!                                                                      *
!     Object:                                                          *
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

use Index_Functions, only: nTri3_Elem
use external_centers, only: iXPolType, nOrd_XF, nXF, nXMolnr, XF, XMolnr
use Symmetry_Info, only: iChBas
use Integral_interfaces, only: int_kernel, int_mem
use Gateway_global, only: PrPrt
use rctfld_module, only: fMax, lMax, lRFCav, Scal14
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: nIrrep
use Basis_Info, only: nBas
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nDens, nCavxyz_, nGrid_, MaxAto
real(kind=wp), intent(in) :: D_Tot(nDens), Grid(3,nGrid_), Cord(3,MaxAto), Z_Nuc(MaxAto)
real(kind=wp), intent(out) :: Ravxyz(nCavxyz_)
real(kind=wp), intent(inout) :: Cavxyz(nCavxyz_), dEF(4,nGrid_), xfEF(4,nGrid_)
integer(kind=iwp) :: iComp, iGrid, iMax, iMltpl, ip, iSym, iSymC, iSymX, iSymXY, iSymXZ, iSymY, iSymYZ, iSymZ, iSyXYZ, iTemp, ix, &
                     ixyz, iy, iz, l_Oper(1), MltLbl, nComp, nh1, nOpr, nOrdOp
real(kind=wp) :: CCoor(3), FactOp(1), fTest, Origin(3), rHrmt, Sig
logical(kind=iwp) :: Save_tmp
character(len=8) :: Label
integer(kind=iwp), allocatable :: ips(:), kOper(:), lOper(:)
real(kind=wp), allocatable :: C_Coor(:,:), Nuc(:)
procedure(int_kernel) :: EFInt
procedure(int_mem) :: EFMem
integer(kind=iwp), external :: IrrFnc
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: lOff, iIrrep, n
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
if (lRFCav) then  !Skip calculation of moments if no cavity
  Origin(:) = Zero
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  ! Compute the multipole moment of the QC system,
  ! both nuclear(1) and electronic(2).

  ! Cavxyz: Multipole moments in cartesian basis
  ! Ravxyz: temporary storage

  Cavxyz(:) = Zero

  ! 1) Compute M(nuc,nl), nuclear multipole moments

  do iMax=0,lMax
    ip = 1+nTri3_Elem(iMax)
    call RFNuc(Origin,Ravxyz(ip),iMax)
  end do

# ifdef _DEBUGPRINT_
  call RecPrt('Nuclear Multipole Moments',' ',Ravxyz,1,nCavxyz_)
# endif

  ! 2) Compute the electronic contribution to the charge distribution.

  ! M(el,nl) =  - Sum(p,q) Dpq <p|M(nl)|q>

  nOpr = 1
  FactOp = One
  l_Oper = 1
  ! Reset array for storage of multipole moment expansion
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
        l_Oper = ior(l_Oper,iTemp)
      end do
    end do
  end do
# ifdef _DEBUGPRINT_
  write(u6,*) '1st order total density'
  lOff = 1
  do iIrrep=0,nIrrep-1
    n = nTri_Elem(nBas(iIrrep))
    write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
    call Triprt(Label,' ',D_Tot(lOff),nBas(iIrrep))
    lOff = lOff+n
  end do
# endif

  call Drv1_RF(FactOp,nOpr,D_tot,nh1,Origin,l_Oper,Cavxyz,lMax)

# ifdef _DEBUGPRINT_
  call RecPrt('Electronic Multipole Moments',' ',Cavxyz,1,nCavxyz_)
# endif

  ! Add nuclear MM expansion to the electronic one

  Cavxyz(:) = Cavxyz(:)+Ravxyz(:)

# ifdef _DEBUGPRINT_
  call RecPrt('Electronic+Nuclear Moments',' ',Cavxyz,1,nCavxyz_)
# endif

  if (allocated(XF)) call XFMoment(lMax,Cavxyz,Ravxyz,nCavxyz_,Origin)

# ifdef _DEBUGPRINT_
  call RecPrt('Total Multipole Moments ',' ',Cavxyz,1,nCavxyz_)
# endif
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over the Langevin grid, compute EF and store in dEF.

! The code here is initially with the loop over the grid outermost.
! We might have to change this later on!

! Be happy, don't worry!

rHrmt = One
nOrdOp = 1
nComp = 3
Sig = -One
iSymX = IrrFnc(1)
iSymY = IrrFnc(2)
iSymZ = IrrFnc(4)
iSymXY = IrrFnc(3)
iSymXZ = IrrFnc(5)
iSymYZ = IrrFnc(6)
iSyXYZ = IrrFnc(7)

call mma_allocate(ips,nComp,Label='ips')
call mma_allocate(lOper,nComp,Label='lOper')
call mma_allocate(kOper,nComp,Label='kOper')
call mma_allocate(Nuc,nComp,Label='Nuc')
call mma_allocate(C_Coor,3,nComp,Label='CCoor')

Save_tmp = PrPrt
PrPrt = .true.

do iGrid=1,nGrid_
  write(Label,'(A,I5)') 'EF ',iGrid
  Ccoor(:) = Grid(:,iGrid)
  iSymC = 1
  if (Ccoor(1) /= Zero) iSymC = ibset(iSymC,iSymX)
  if (Ccoor(2) /= Zero) iSymC = ibset(iSymC,iSymY)
  if (Ccoor(3) /= Zero) iSymC = ibset(iSymC,iSymZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero)) iSymC = ibset(iSymC,iSymXY)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ibset(iSymC,iSymXZ)
  if ((Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ibset(iSymC,iSymYZ)
  if ((Ccoor(1) /= Zero) .and. (Ccoor(2) /= Zero) .and. (Ccoor(3) /= Zero)) iSymC = ibset(iSymC,iSyXYZ)

  iComp = 0
  do ix=nOrdOp,0,-1
    do iy=nOrdOp-ix,0,-1
      iComp = iComp+1
      iz = nOrdOp-ix-iy
      ixyz = 0
      if (mod(ix,2) /= 0) ixyz = ibset(ixyz,0)
      if (mod(iy,2) /= 0) ixyz = ibset(ixyz,1)
      if (mod(iz,2) /= 0) ixyz = ibset(ixyz,2)
      iSym = 2**IrrFnc(ixyz)
      if (Ccoor(iComp) /= Zero) iSym = ibset(iSym,0)
      lOper(iComp) = MltLbl(iSymC,iSym)
      kOper(iComp) = iChBas(iComp+1)
      C_Coor(:,iComp) = Ccoor(:)
    end do
  end do

  call EFNuc(C_Coor,Z_Nuc,Cord,MaxAto,Nuc,nOrdOp)
  call OneEl_Property(EFInt,EFMem,Label,ips,lOper,nComp,C_Coor,nOrdOp,Nuc,rHrmt,kOper,D_Tot,nDens,dEF(:,iGrid),Sig)

  ! Field contribution from XF
  call EFXF(C_Coor,XF,nXF,nOrd_XF,iXPolType,xfEF(:,iGrid),XMolnr,nXMolnr,iGrid,scal14)

end do

! Add XF contribution to the total field
dEF(:,:) = dEF(:,:)+xfEF(:,:)

PrPrt = Save_tmp

call mma_deallocate(C_Coor)
call mma_deallocate(Nuc)
call mma_deallocate(kOper)
call mma_deallocate(lOper)
call mma_deallocate(ips)
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over grid points and compute square norm of EF. Store the
! largest of the latter in fmax (used in solver).

fmax = Zero
do iGrid=1,nGrid_
  ftest = dEF(1,iGrid)**2+dEF(2,iGrid)**2+dEF(3,iGrid)**2
  dEF(4,iGrid) = ftest
  fmax = max(fmax,ftest)
end do ! iGrid

!call RecPrt('eperm: dEF ',' ',dEF,4,nGrid_)

return

end subroutine eperm
