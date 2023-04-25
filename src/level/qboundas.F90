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
subroutine QBOUNDas(KV,JROT,E,EO,VMX,DSOC,VBZ,SDRDY,RVB,YMIN,YH,GB,GI,SB,SI,NPP,ITP2,ITP3,IWR,IQTST,BFCT,IT)
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

integer I, II, IQTST, IT, ITP2, ITP3, IWR, J, JROT, KV, NPP
real*8 A1, A2, A13, A23, BFCT, C1A, C2A, DSOC, E, EO, FBA, FIA, FJ, GB, GI, GBA, GIA, YH, RH, YMIN, YMINN, SB, SI, SL, VBZ(NPP), &
       SDRDY(NPP), RVB(NPP), VMX, VMXPR
data C1A/0.355028053887817d0/,C2A/0.258819403792807d0/

IQTST = 1
YMINN = YMIN-YH
write(6,*) 'YMINN=',YMINN ! Make sure YMINN is referenced.
! Start by searching for third turning point.
J = NPP-1
if (VBZ(J) > E) go to 30
do I=NPP-2,1,-1
  J = J-1
  if (VBZ(J) > E) go to 10
end do
IQTST = -9
write(6,602) JROT,EO
return
! ITP3 is the first mesh point outside classically forbidden region
10 continue
ITP3 = J+1
! Check that there is a classically allowed region inside this point
! and determine height of barrier maximum.
II = J
VMX = DSOC
do I=2,J
  II = II-1
  if (VBZ(II) <= E) go to 20
  if (VBZ(II) > VMX) VMX = VBZ(II)
end do
! Energy too high (or too low): find only one turning point.
VMXPR = VMX/BFCT
if (IWR /= 0) write(6,604) JROT,EO,VMXPR/BFCT,RVB(J)
IQTST = -1
return
! ITP2 is first mesh point inside forbidden region on left of barrier
20 continue
ITP2 = II+1
! Now ... continue to set up r3(E) boundary condition ...
RH = RVB(ITP3)-RVB(ITP3-1)
GB = (VBZ(ITP3)-E)*(RH/YH)**2
GI = (VBZ(ITP3-1)-E)*(RH/YH)**2
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING
!write(6,*) 'ITP3=',ITP3
!write(6,*) 'VBZ(ITP3-1)=',VBZ(ITP3-1)
!write(6,*) 'E=',E
!write(6,*) 'RH=',RH
!write(6,*) 'YH=',YH
FJ = GI/(GI-GB)
write(6,*) FJ ! make sure it's "referenced"
! Treat quasibound levels as bound using outer boundary condition
! of Airy function at third turning point ... as discussed by
! R.J.Le Roy and R.B.Bernstein  in  J.Chem.Phys. 54,5114(1971).
! Uses series expansions of Abramowitz & Stegun Eq.(10.4.3)
SL = (GI-GB)**(1.d0/3.d0)/RH
A1 = GI/(SL*RH)**2
A2 = GB/(SL*RH)**2
A13 = A1*A1*A1
A23 = A2*A2*A2
FIA = 1.d0+A13*(A13*(A13+72.d0)+2160.d0)/12960.d0
GIA = A1+A1*A13*(A13*(A13+90.d0)+3780.d0)/45360.d0
FBA = 1.d0+A23*(A23*(A23+72.d0)+2160.d0)/12960.d0
GBA = A2+A2*A23*(A23*(A23+90.d0)+3780.d0)/45360.d0
! Airy function  Bi(X)  at points straddling 3-rd turning point
SI = (C1A*FIA+C2A*GIA)/SDRDY(ITP3-1)
SB = (C1A*FBA+C2A*GBA)/SDRDY(ITP3)
GI = VBZ(ITP3-1)-E
GB = VBZ(ITP3)-E
if (SB >= SI) then
  ! In case of big error - switch to node at ITP3
  SB = 0.d0
  SI = 1.d0
  if (IWR /= 0) write(6,606) KV,JROT,EO,IT
end if
return

! If 3-rd turning point beyond range start with WKB wave function
! at end of range.
30 continue
if (IWR /= 0) write(6,608) JROT,EO
ITP3 = NPP-1
IQTST = 0
VMX = VBZ(ITP3)
II = ITP3
!... and determine barrier maximum ....
do I=2,ITP3
  II = II-1
  VMXPR = VBZ(II)
  if (VMXPR < VMX) go to 40
  VMX = VMXPR
end do
if (IWR /= 0) write(6,610)
IQTST = -9

40 continue
return

602 format(' *** QBOUND fails for   E(J=',i3,')=',f9.3,'  Find no turning point')
604 format(' For J=',I3,'  ETRY=',F11.4,' > VMAX=',F11.4,'  find onee turn point:  R=',F6.2)
606 format(' *** CAUTION ***  v=',I3,'   J=',I3,'   E=',1PD13.6,'   IT=',I2/5x, &
           'Airy initialization unstable so place node just past  R(3-rd)')
608 format(' *** For  J=',I3,'  E=',F9.2,'  R(3-rd) > RMAX  & E < V(N)  so try WKB B.C. @ RMAX')
610 format(" **** QBOUND doesn't work ... no classically allowed region accessible at this energy.")

end subroutine QBOUNDas
