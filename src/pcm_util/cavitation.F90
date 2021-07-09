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

subroutine Cavitation(DoDeriv,ToAng,nAt,nS,nTs,GCavP,VMol,TAbs,RSolv,Sphere,Tessera,iSphe)
! Compute cavitation energy given the cavity.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: DoDeriv
real(kind=wp), intent(in) :: ToAng, VMol, TAbs, RSolv
integer(kind=iwp), intent(in) :: nAt, nS, nTs, iSphe(nTs)
real(kind=wp), intent(out) :: GCavP
real(kind=wp), intent(inout) :: Sphere(4,nS), Tessera(4,nTs)
integer(kind=iwp) :: iS, iTs
real(kind=wp) :: ToAU
real(kind=wp), allocatable :: CavS(:), dCav(:,:), dEA(:,:,:), EA(:)

call mma_allocate(CavS,nS,label='CavSph')
call mma_allocate(dCav,3,nAt,label='dCav')
call mma_allocate(EA,nS,label='ExpArea')
call mma_allocate(dEA,3,nAt,nS,label='dExpArea')
CavS(:) = Zero
dCav(:,:) = Zero
EA(:) = Zero
dEA(:,:,:) = Zero

! Put some quantities in angstrom
call DScal_(nS,ToAng,Sphere(4,1),4)
call DScal_(nTs,ToAng*ToAng,Tessera(4,1),4)

! Compute the exposed area for each sphere
do iTs=1,nTs
  iS = iSphe(iTs)
  EA(iS) = EA(iS)+Tessera(4,iTs)
end do

call Cavitation_(nAt,nS,VMol,TAbs,RSolv,GCavP,CavS,dCav,Sphere,EA,dEA,DoDeriv)

! Revert to a.u.
ToAU = One/ToAng
call DScal_(nS,ToAU,Sphere(4,1),4)
call DScal_(nTs,ToAU*ToAU,Tessera(4,1),4)

call mma_deallocate(CavS)
call mma_deallocate(dCav)
call mma_deallocate(EA)
call mma_deallocate(dEA)

return

end subroutine Cavitation
!====
subroutine Cavitation_(NAtoms,NSph,VMol,TAbs,RSolv,GCavP,PCvSph,dCav,Sphere,Ae,dAe,DoDeriv)
! Compute cavitation energy given the cavity.

use Constants, only: Zero, One, Two, Three, Pi, rNAVO, Rgas, cal_to_J
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NAtoms, NSph
real(kind=wp), intent(in) :: VMol, TAbs, RSolv, Sphere(4,NSph), Ae(NSph), dAe(3,NAtoms,NSph)
real(kind=wp), intent(out) :: GCavP, PCvSph(NSph)
real(kind=wp), intent(inout) :: dCav(3,NAtoms)
logical(kind=iwp), intent(in) :: DoDeriv
integer(kind=iwp) :: IAtom, ISph, Ixyz
real(kind=wp) :: AExpsd, ATotal, DENum, Fact, FPI, Frac, G, R, RSph, RT, TPI, YM, YP
real(kind=wp), parameter :: Avgdr = rNAVO*1.0e-24_wp, GC = Rgas/cal_to_J

TPI = Two*PI
FPI = Two*TPI

! Cavitation Energy:  R.A. Pierotti, Chem.Rev. 76,717,(1976).
!
! The cavitation energy is now computed for each sphere of the
! cavity. The values are scaled for the exposed surface of each
! sphere and summed irrespective of the number of cavities.
! The derivative of the cavitation energy wrt nuclear coordinates
! is also computed. To this end the Pierotti cavitation energy for
! each sphere is stored in the array PCvSph(ISph).

DENum = Avgdr/VMol
RT = GC*TAbs*0.001_wp
YP = DENum*FPI*RSolv**3/Three
YM = YP/(One-YP)

! Loop on spheres.

GCavP = Zero
do ISph=1,NSph
  RSph = Sphere(4,ISph)
  R = RSph/RSolv

  ! Pierotti's cavitation free energy for the single sphere.

  G = -log(One-YP)+4.5_wp*YM**2*R**2+Three*YM*R*(R+One)
  G = RT*G
  PCvSph(ISph) = G
  AExpsd = Ae(ISph)
  ATotal = FPI*RSph*RSph
  Frac = AExpsd/ATotal
  GCavP = GCavP+G*Frac
end do

! Computation of derivatives wrt nuclear coordinates according to
! the modification of Claverie to the Pierotti approach.

! Loop on atoms, coordinates, tesserae and then spheres.

if (DoDeriv) then
  do ISph=1,NSph
    RSph = Sphere(4,ISph)
    Fact = PCvSph(ISph)/(FPI*RSph*RSph)
    do IAtom=1,NAtoms
      do Ixyz=1,3
        dCav(Ixyz,IAtom) = dCav(Ixyz,IAtom)+Fact*dAe(Ixyz,IAtom,ISph)
      end do
    end do
  end do
end if

return

end subroutine Cavitation_
