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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine QBOUNDas(KV,JROT,E,EO,VMX,DSOC,VBZ,SDRDY,RVB,YH,GB,GI,SB,SI,NPP,ITP2,ITP3,IWR,IQTST,BFCT,IT)
!***********************************************************************
!** Subroutine to initialize quasibound level wave function as Airy
!  function at third turning point (if possible). For the theory see
!  J.Chem.Phys. 54, 5114 (1971),  J.Chem.Phys. 69, 3622-31 (1978)
!----------------------------------------------------------------------
!** IQTST  is error flag. *** If (IQTST < 0) initialization fails
!  so eigenvalue calculation aborts *** (IQTST > 0) for successful
!  Airy function initialization. *** (IQTST=0) if Airy function
!  initialization prevented because 3-rd turning point beyond
!  range, so that WKB initialization is used.
!----------------------------------------------------------------------

use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: KV, JROT, NPP, IWR, IT
real(kind=wp), intent(in) :: E, EO, DSOC, VBZ(NPP), SDRDY(NPP), RVB(NPP), YH, BFCT
integer(kind=iwp), intent(out) :: ITP2, ITP3, IQTST
real(kind=wp), intent(out) :: VMX, GB, GI, SB, SI
integer(kind=iwp) :: I, II, J
real(kind=wp) :: A1, A13, A2, A23, FBA, FIA, GBA, GIA, RH, SL, VMXPR
logical(kind=iwp) :: Found
#include "compiler_features.h"
#define _C1A_ One/(Three**(Two/Three)*gamma(Two/Three))
#define _C2A_ One/(Three**(One/Three)*gamma(One/Three))
#ifdef INTRINSIC_INITIALIZATION
real(kind=wp), parameter :: C1A = _C1A_, C2A = _C2A_
#else
real(kind=wp) :: C1A, C2A
C1A = _C1A_
C2A = _C2A_
#endif

IQTST = 1
! Start by searching for third turning point.
J = NPP-1
if (VBZ(J) > E) then
  ! If 3-rd turning point beyond range start with WKB wave function
  ! at end of range.
  if (IWR /= 0) write(u6,608) JROT,EO
  ITP3 = NPP-1
  IQTST = 0
  VMX = VBZ(ITP3)
  II = ITP3
  !... and determine barrier maximum ....
  Found = .false.
  do I=2,ITP3
    II = II-1
    VMXPR = VBZ(II)
    if (VMXPR < VMX) then
      Found = .true.
      exit
    end if
    VMX = VMXPR
  end do
  if (.not. Found) then
    if (IWR /= 0) write(u6,610)
    IQTST = -9
  end if
  return
end if
Found = .false.
do I=NPP-2,1,-1
  J = J-1
  if (VBZ(J) > E) then
    Found = .true.
    exit
  end if
end do
if (.not. Found) then
  IQTST = -9
  write(u6,602) JROT,EO
  return
end if
! ITP3 is the first mesh point outside classically forbidden region
ITP3 = J+1
! Check that there is a classically allowed region inside this point
! and determine height of barrier maximum.
II = J
VMX = DSOC
Found = .false.
do I=2,J
  II = II-1
  if (VBZ(II) <= E) then
    Found = .true.
    exit
  end if
  if (VBZ(II) > VMX) VMX = VBZ(II)
end do
if (.not. Found) then
  ! Energy too high (or too low): find only one turning point.
  VMXPR = VMX/BFCT
  if (IWR /= 0) write(u6,604) JROT,EO,VMXPR/BFCT,RVB(J)
  IQTST = -1
  return
end if
! ITP2 is first mesh point inside forbidden region on left of barrier
ITP2 = II+1
! Now ... continue to set up r3(E) boundary condition ...
RH = RVB(ITP3)-RVB(ITP3-1)
GB = (VBZ(ITP3)-E)*(RH/YH)**2
GI = (VBZ(ITP3-1)-E)*(RH/YH)**2
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING
!write(u6,*) 'ITP3=',ITP3
!write(u6,*) 'VBZ(ITP3-1)=',VBZ(ITP3-1)
!write(u6,*) 'E=',E
!write(u6,*) 'RH=',RH
!write(u6,*) 'YH=',YH
! Treat quasibound levels as bound using outer boundary condition
! of Airy function at third turning point ... as discussed by
! R.J.Le Roy and R.B.Bernstein  in  J.Chem.Phys. 54,5114(1971).
! Uses series expansions of Abramowitz & Stegun Eq.(10.4.3)
SL = (GI-GB)**(One/Three)/RH
A1 = GI/(SL*RH)**2
A2 = GB/(SL*RH)**2
A13 = A1*A1*A1
A23 = A2*A2*A2
FIA = One+A13*(A13*(A13+72.0_wp)+2160.0_wp)/12960.0_wp
GIA = A1+A1*A13*(A13*(A13+90.0_wp)+3780.0_wp)/45360.0_wp
FBA = One+A23*(A23*(A23+72.0_wp)+2160.0_wp)/12960.0_wp
GBA = A2+A2*A23*(A23*(A23+90.0_wp)+3780.0_wp)/45360.0_wp
! Airy function  Bi(X)  at points straddling 3-rd turning point
SI = (C1A*FIA+C2A*GIA)/SDRDY(ITP3-1)
SB = (C1A*FBA+C2A*GBA)/SDRDY(ITP3)
GI = VBZ(ITP3-1)-E
GB = VBZ(ITP3)-E
if (SB >= SI) then
  ! In case of big error - switch to node at ITP3
  SB = Zero
  SI = One
  if (IWR /= 0) write(u6,606) KV,JROT,EO,IT
end if

return

602 format(' *** QBOUND fails for   E(J=',i3,')=',f9.3,'  Find no turning point')
604 format(' For J=',I3,'  ETRY=',F11.4,' > VMAX=',F11.4,'  find onee turn point:  R=',F6.2)
606 format(' *** CAUTION ***  v=',I3,'   J=',I3,'   E=',ES13.6,'   IT=',I2/5x, &
           'Airy initialization unstable so place node just past  R(3-rd)')
608 format(' *** For  J=',I3,'  E=',F9.2,'  R(3-rd) > RMAX  & E < V(N)  so try WKB B.C. @ RMAX')
610 format(" **** QBOUND doesn't work ... no classically allowed region accessible at this energy.")

end subroutine QBOUNDas
