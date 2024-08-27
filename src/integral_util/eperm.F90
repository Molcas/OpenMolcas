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

use external_centers, only: iXPolType, nOrd_XF, nXF, nXMolnr, XF, XMolnr
use Symmetry_Info, only: iChBas
#ifdef _DEBUGPRINT_
use Symmetry_Info, only: nIrrep
use Basis_Info, only: nBas
#endif
use Integral_interfaces, only: int_kernel, int_mem
use Gateway_global, only: PrPrt
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use rctfld_module, only: lRFCav, lMax, nGrid, fMax, Scal14

implicit none
integer nDens, nCavxyz_, nGrid_, MaxAto
real*8 D_Tot(nDens), Ravxyz(nCavxyz_), Cavxyz(nCavxyz_), dEF(4,nGrid_), Grid(3,nGrid_), Cord(3,MaxAto), Z_Nuc(MaxAto), &
       xfEF(4,nGrid_)
real*8 Origin(3), CCoor(3)
procedure(int_kernel) :: EFInt
procedure(int_mem) :: EFMem
logical Save_tmp
character*8 Label
real*8 FactOp(1)
integer l_Oper(1)
integer, allocatable :: ips(:), lOper(:), kOper(:)
real*8, allocatable :: C_Coor(:,:), Nuc(:)
integer ixyz, iOff
integer iMax, ip, iMltpl, ix, iy, iz, iSymX, iSymY, iSymZ, iTemp, nComp, iSymXY, iSymXZ, iSymYZ, iSyXYZ, iComp, iSym, nh1, MltLbl, &
        nOpr, nOrdOp, iGrid, iSymC
integer, external :: IrrFnc
real*8 rHrmt, Sig, fTest
#ifdef _DEBUGPRINT_
integer lOff, iIrrep, n
#endif
! Statement Function
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. lRFCav) goto 99  !Skip calculation of moments if no cavity
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
  ip = 1+iOff(iMax)
  call RFNuc(Origin,Ravxyz(ip),iMax)
end do

#ifdef _DEBUGPRINT_
call RecPrt('Nuclear Multipole Moments',' ',Ravxyz,1,nCavxyz_)
#endif

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
      l_Oper = ior(l_Oper,iTemp)
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(6,*) '1st order total density'
lOff = 1
do iIrrep=0,nIrrep-1
  n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
  write(Label,'(A,I1)') 'Diagonal Symmetry Block ',iIrrep+1
  call Triprt(Label,' ',D_Tot(lOff),nBas(iIrrep))
  lOff = lOff+n
end do
#endif

call Drv1_RF(FactOp,nOpr,D_tot,nh1,Origin,l_Oper,Cavxyz,lMax)

#ifdef _DEBUGPRINT_
call RecPrt('Electronic Multipole Moments',' ',Cavxyz,1,nCavxyz_)
#endif

! Add nuclear MM expansion to the electronic one

call DaXpY_(nCavxyz_,One,Ravxyz,1,Cavxyz,1)

#ifdef _DEBUGPRINT_
call RecPrt('Electronic+Nuclear Moments',' ',Cavxyz,1,nCavxyz_)
#endif

if (allocated(XF)) call XFMoment(lMax,Cavxyz,Ravxyz,nCavxyz_,Origin)

#ifdef _DEBUGPRINT_
call RecPrt('Total Multipole Moments ',' ',Cavxyz,1,nCavxyz_)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over the Langevin grid, compute EF and store in dEF.

! The code here is initially with the loop over the grid outermost.
! We might have to change this later on!

! Be happy, don't worry!

99 continue
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
call mma_allocate(Nuc,nComp,Label='Nuc')
call mma_allocate(C_Coor,3,nComp,Label='CCoor')

Save_tmp = PrPrt
PrPrt = .true.

do iGrid=1,nGrid_
  write(Label,'(A,I5)') 'EF ',iGrid
  call dcopy_(3,Grid(1,iGrid),1,Ccoor,1)
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

  call EFNuc(C_Coor,Z_Nuc,Cord,MaxAto,Nuc,nOrdOp)
  call OneEl_Property(EFInt,EFMem,Label,ips,lOper,nComp,C_Coor,nOrdOp,Nuc,rHrmt,kOper,D_Tot,nDens,dEF(1,iGrid),Sig)

  ! Field contribution from XF
  call EFXF(C_Coor,XF,nXF,nOrd_XF,iXPolType,xfEF(1,iGrid),XMolnr,nXMolnr,iGrid,scal14)

end do

! Add XF contribution to the total field
call DaXpY_(4*nGrid,One,xfEF,1,dEF,1)

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
