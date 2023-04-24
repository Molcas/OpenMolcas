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
!
!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
      SUBROUTINE QBOUNDas(KV,JROT,E,EO,VMX,DSOC,VBZ,SDRDY,RVB,YMIN,YH,  &
     &  GB,GI,SB,SI,NPP,ITP2,ITP3,IWR,IQTST,BFCT,IT)
!***********************************************************************
!** Subroutine to initialize quasibound level wave function as Airy
!  function at third turning point (if possible). For the theory see
!  J.Chem.Phys. 54, 5114 (1971),  J.Chem.Phys. 69, 3622-31 (1978)
!----------------------------------------------------------------------
!** IQTST  is error flag. *** If (IQTST.lt.0) initialization fails
!  so eigenvalue calculation aborts *** (IQTST.gt.0) for successful
!  Airy function initialization. *** (IQTST=0) if Airy function
!  initialization prevented because 3-rd turning point beyond
!  range, so that WKB initialization is used.
!----------------------------------------------------------------------
      INTEGER I,II,IQTST,IT,ITP2,ITP3,IWR,J,JROT,KV,NPP
      REAL*8  A1,A2,A13,A23,BFCT,C1A,C2A,DSOC,E,EO,FBA,FIA,FJ,GB,GI,    &
     &  GBA,GIA,YH,RH,YMIN,YMINN,SB,SI,SL,VBZ(NPP),SDRDY(NPP),RVB(NPP), &
     &  VMX,VMXPR
      DATA C1A/0.355028053887817D0/,C2A/0.258819403792807D0/
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IQTST=1
      YMINN= YMIN- YH
      WRITE(6,*) 'YMINN=',YMINN ! Make sure YMINN is referenced.
!** Start by searching for third turning point.
      J=NPP-1
      IF(VBZ(J).GT.E) GO TO 30
      DO  I=NPP-2,1,-1
          J=J-1
          IF(VBZ(J).GT.E) GO TO 10
          ENDDO
      IQTST= -9
      WRITE(6,602) JROT,EO
      RETURN
!** ITP3 is the first mesh point outside classically forbidden region
   10 ITP3= J+1
!** Check that there is a classically allowed region inside this point
!  and determine height of barrier maximum.
      II=J
      VMX=DSOC
      DO  I=2,J
          II=II-1
          IF(VBZ(II).LE.E) GO TO 20
          IF(VBZ(II).GT.VMX) VMX= VBZ(II)
          ENDDO
!** Energy too high (or too low): find only one turning point.
      VMXPR= VMX/BFCT
      IF(IWR.NE.0) WRITE(6,604) JROT,EO,VMXPR/BFCT,RVB(J)
      IQTST=-1
      RETURN
!** ITP2 is first mesh point inside forbidden region on left of barrier
   20 ITP2= II+1
!** Now ... continue to set up r3(E) boundary condition ...
      RH= RVB(ITP3)- RVB(ITP3-1)
      GB= (VBZ(ITP3) - E)*(RH/YH)**2
      GI= (VBZ(ITP3-1) - E)*(RH/YH)**2
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING
!     WRITE(6,*) 'ITP3=',ITP3
!     WRITE(6,*) 'VBZ(ITP3-1)=',VBZ(ITP3-1)
!     WRITE(6,*) 'E=',E
!     WRITE(6,*) 'RH=',RH
!     WRITE(6,*) 'YH=',YH
      FJ= GI/(GI-GB)
      WRITE(6,*) FJ ! make sure it's "referenced"
!** Treat quasibound levels as bound using outer boundary condition
!  of Airy function at third turning point ... as discussed by
!  R.J.Le Roy and R.B.Bernstein  in  J.Chem.Phys. 54,5114(1971).
!  Uses series expansions of Abramowitz & Stegun Eq.(10.4.3)
      SL= (GI-GB)**(1.d0/3.d0)/RH
      A1= GI/(SL*RH)**2
      A2= GB/(SL*RH)**2
      A13= A1*A1*A1
      A23= A2*A2*A2
      FIA= 1.d0+ A13*(A13*(A13+72.D0)+2160.D0)/12960.D0
      GIA= A1+A1*A13*(A13*(A13+90.D0)+3780.D0)/45360.D0
      FBA= 1.d0+ A23*(A23*(A23+72.D0)+2160.D0)/12960.D0
      GBA= A2+A2*A23*(A23*(A23+90.D0)+3780.D0)/45360.D0
!** Airy function  Bi(X)  at points straddling 3-rd turning point
      SI= (C1A*FIA+C2A*GIA)/SDRDY(ITP3-1)
      SB= (C1A*FBA+C2A*GBA)/SDRDY(ITP3)
      GI= VBZ(ITP3-1) - E
      GB= VBZ(ITP3) - E
      IF(SB.GE.SI) THEN
!** In case of big error - switch to node at ITP3
          SB= 0.d0
          SI= 1.d0
          IF(IWR.NE.0) WRITE(6,606) KV,JROT,EO,IT
          ENDIF
      RETURN
!
!** If 3-rd turning point beyond range start with WKB wave function
!  at end of range.
   30 IF(IWR.NE.0) WRITE(6,608) JROT,EO
      ITP3= NPP-1
      IQTST= 0
      VMX= VBZ(ITP3)
      II= ITP3
!... and determine barrier maximum ....
      DO  I= 2,ITP3
          II= II-1
          VMXPR= VBZ(II)
          IF(VMXPR.LT.VMX) GO TO 40
          VMX= VMXPR
          ENDDO
      IF(IWR.NE.0) WRITE(6,610)
      IQTST= -9
   40 RETURN
  602 FORMAT(' *** QBOUND fails for   E(J=',i3,')=',f9.3,'  Find no turn&
     &ing point')
  604 FORMAT(' For J=',I3,'  ETRY=',F11.4,' > VMAX=',F11.4,             &
     &  '  find onee turn point:  R=',F6.2)
  606 FORMAT(' *** CAUTION ***  v=',I3,'   J=',I3,'   E=',1PD13.6,      &
     & '   IT=',I2/5x,'Airy initialization unstable so place node just p&
     &ast  R(3-rd)' )
  608 FORMAT(' *** For  J=',I3,'  E=',F9.2,                             &
     &  '  R(3-rd) > RMAX  & E < V(N)  so try WKB B.C. @ RMAX')
  610 FORMAT(" **** QBOUND doesn't work ... no classically allowed regio&
     &n accessible at this energy.")
      END
