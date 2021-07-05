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

implicit real*8(A-H,O-Z)
logical DoDeriv
#include "WrkSpc.fh"
dimension Sphere(4,*), Tessera(4,*), iSphe(*)

call GetMem('CavSph','Allo','Real',ip_CavS,nS)
call GetMem('dCav','Allo','Real',ip_dCav,3*nAt)
call GetMem('ExpArea','Allo','Real',ip_EA,nS)
call GetMem('dExpArea','Allo','Real',ip_dEA,3*nAt*nS)
call FZero(Work(ip_CavS),nS)
call FZero(Work(ip_dCav),3*nAt)
call FZero(Work(ip_EA),nS)
call FZero(Work(ip_dEA),3*nAt*nS)

! Put some quantity in Angstrom
call DScal_(nS,ToAng,Sphere(4,1),4)
call DScal_(nTs,ToAng*ToAng,Tessera(4,1),4)

! Compute the exposed area for each sphere
do iTs=1,nTs
  iS = iSphe(iTs)
  Work(ip_EA-1+iS) = Work(ip_EA-1+iS)+Tessera(4,iTs)
end do

call Cavitation_(nAt,nS,VMol,TAbs,RSolv,GCavP,Work(ip_CavS),Work(ip_dCav),Sphere,Work(ip_EA),Work(ip_dEA),DoDeriv)

! Revert to a.u.
ToAU = 1.0d0/ToAng
call DScal_(nS,ToAU,Sphere(4,1),4)
call DScal_(nTs,ToAU*ToAU,Tessera(4,1),4)

call GetMem('dExpArea','Free','Real',ip_dEA,3*nAt*nS)
call GetMem('ExpArea','Free','Real',ip_EA,nS)
call GetMem('dCav','Free','Real',ip_dCav,3*nAt)
call GetMem('CavSph','Free','Real',ip_CavS,nS)

return

end subroutine Cavitation
!====
subroutine Cavitation_(NAtoms,NSph,VMol,TAbs,RSolv,GCavP,PCvSph,dCav,Sphere,Ae,dAe,DoDeriv)
! Compute cavitation energy given the cavity.

implicit real*8(A-H,O-Z)
logical DoDeriv
real*8 Sphere(4,*), Ae(*), dAe(3,NAtoms,*), PCvSph(*), dCav(3,*)
save Zero, One, Two, Three, Four, F4Pt5, F1000, GC, Avgdr
data Zero,One,Two,Three,Four,F4Pt5,F1000/0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,4.5d0,1.0d3/
data GC,Avgdr/1.9865d0,0.60228d0/

PI = Four*atan(One)
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
RT = GC*TAbs/F1000
YP = DENum*FPI*RSolv**3/Three
YM = YP/(One-YP)

! Loop on spheres.

GCavP = Zero
do ISph=1,NSph
  RSph = Sphere(4,ISph)
  R = RSph/RSolv

  ! Pierotti's cavitation free energy for the single sphere.

  G = -log(One-YP)+F4Pt5*YM**2*R**2+Three*YM*R*(R+One)
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
