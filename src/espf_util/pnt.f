************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Christophe Chipot                                      *
*               Janos G. Angyan                                        *
************************************************************************
      Subroutine PNT(IOut,nAtoms,CO,IRMax,DeltaR,IAn,nPts,Grid,
     &               IsMM,Process)
      Implicit Real*8(A-H,O-Z)
C
c     THIS SUBROUTINE CREATES A GRID OF POINTS AROUND A MOLECULE.
C
c     THE POINTS ARE INITIALLY SELECTED IN A CUBE, THEN RE-SCALED TO THE
c     SIZE OF THE MOLECULE.
C
c     THE POINTS LYING WITHIN THE VAN DER WAALS ENVELOPE OF THE
c     MOLECULE ARE EXCLUDED, BECAUSE OF THE PENETRATION EFFECTS.
c     THIS PROCEDURE IS BASED ON THE SCALED VDW RADII.
c     THE SCALING FACTOR : VDWFact
C
c     THE POINTS LYING OUTSIDE THE EXCLUSION RADIUS ARE EXCLUDED TOO.
c     THE EXCLUSION RADIUS  :  r_max
c     THE GRID STEP         :  delta_r
C
C     This subroutine comes from the Grid_3.1.2 program (C.CHIPOT, J.G.ANGYAN,
C     Universite' Henri Poincare' - Nancy 1, France)
C
      Logical Retain1,Retain2,Process
      Integer IsMM(nAtoms)
      Dimension CO(3,nAtoms),IAn(nAtoms),Grid(3,*),RVDW(0:54)
      Save Zero,Two,Huge,RVDW
      Data Zero/0.0D0/,Two/2.0D0/,Huge/1.0D8/,
     &     RVDW/2.00D0,1.20D0,1.40D0,
     1          1.82D0,0.00D0,0.00D0,1.70D0,1.55D0,1.50D0,1.47D0,
     2          1.54D0,2.27D0,1.73D0,2.30D0,2.10D0,1.80D0,1.80D0,
     3          1.75D0,1.88D0,2.75D0,0.00D0,0.00D0,0.00D0,1.22D0,
     4          0.00D0,0.00D0,0.00D0,0.00D0,1.63D0,1.40D0,1.39D0,
     5          1.87D0,0.00D0,1.85D0,1.90D0,1.85D0,2.02D0,0.00D0,
     6          0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.0D00,
     7          0.00D0,1.63D0,1.72D0,1.58D0,1.93D0,2.17D0,0.00D0,
     8          2.06D0,1.98D0,2.16D0/
C
c     Criteria ?
C
      iPL = iPL_espf()-1
      RMax=dble(IRMax)
      If (Process .and. iPL.ge.3) then
        Write(IOut,'(A,I2,A)') ' Max : ',IRMax,' van der Waals radii'
        Write(IOut,'(A,F4.2,A)') ' ... with ',DeltaR,' angstroms'//
     &                           ' between grid points.'
      End If
C
c     The VDW scaling factor is currently hard-coded. Too be changed.
C
      VDWFact = 1.0d0
C
C
c     SEARCH FOR THE EXTREMA OF THE MOLECULAR GEOMETRY WITHIN A CUBE.
C
      xmax = - Huge
      xmin =   Huge
      ymax = - Huge
      ymin =   Huge
      zmax = - Huge
      zmin =   Huge
C
      RVDWMX = Zero
      Do 10 iAtom = 1, nAtoms
         If(IsMM(iAtom).ne.0) Then
           If (Process .and. iPL.ge.3) then
              Write(6,11) iAtom
11         Format (' MM atom',I3,' is ignored in the grid construction')
           End If
           Goto 10
         End If
         xmax   = Max(xmax,CO(1,iAtom))
         xmin   = Min(xmin,CO(1,iAtom))
         ymax   = Max(ymax,CO(2,iAtom))
         ymin   = Min(ymin,CO(2,iAtom))
         zmax   = Max(zmax,CO(3,iAtom))
         zmin   = Min(zmin,CO(3,iAtom))
         RVDWMX = Max(RVDWMX,RVDW(IAn(iAtom)))
10    Continue
      rSpher = RMax*RVDWMX
C
c     DETERMINE THE MINIMUM PARALLELEPIPED DIMENSIONS REQUIRED TO
c     CONTAIN THE MOLECULE, INCLUDING A MAXIMUM SELECTION RADIUS
c     ON BOTH SIDES.
C
      xrange = xmax - xmin + Two*rSpher
      yrange = ymax - ymin + Two*rSpher
      zrange = zmax - zmin + Two*rSpher
C
      nbptx = Int(xrange/DeltaR)
      nbpty = Int(yrange/DeltaR)
      nbptz = Int(zrange/DeltaR)
C
C
c     DETERMINE THE MAXIMUM NUMBER OF POINTS FOR THE FITTING PROCEDURE.
C
      maxpt = nbptx * nbpty * nbptz
C
C
c     PRINT OUT INTERMEDIATE RESULTS, IF SO REQUESTED
C
      If (iPL.gt.2) then
         Write(IOut,'(" EXTREMA OF THE MOLECULAR GEOMETRY : ")')
         Write(IOut,'(a)')
         Write(IOut,'(" X_min = ",f8.4,4x,"X_max = ",f8.4)')
     &        xmin, xmax
         Write(IOut,'(" Y_min = ",f8.4,4x,"Y_max = ",f8.4)')
     &        ymin, ymax
         Write(IOut,'(" Z_min = ",f8.4,4x,"Z_max = ",f8.4)')
     &        zmin, zmax
         Write(IOut,'(a)')
         Write(IOut,'(" RE-SCALED PARALLELEPIPED SIZE : ")')
         Write(IOut,'(a)')
         Write(IOut,'(" X = ",f8.4)') xrange
         Write(IOut,'(" Y = ",f8.4)') yrange
         Write(IOut,'(" Z = ",f8.4)') zrange
         Write(IOut,'(a)')
         Write(IOut,'(" NUMBER OF POINTS PER DIRECTION : ")')
         Write(IOut,'(a)')
         Write(IOut,'(" Nb_x = ",i5)') nbptx
         Write(IOut,'(" Nb_y = ",i5)') nbpty
         Write(IOut,'(" Nb_z = ",i5)') nbptz
         Write(IOut,'(a)')
         Write(IOut,'(" MAXIMUM POSSIBLE NUMBER OF POINTS : ")')
         Write(IOut,'(a)')
         Write(IOut,'(" Nb_max = ",i10)') maxpt
         Write(IOut,'(a)')
      End If
C
C
c     LOOP OVER ALL POSSIBLE POINTS IN THE PARALLELEPIPED.
C
      nPts = 0
C
      Do 20 Ipx = 0, nbptx
         pointa = xmin - rSpher + Dble(Ipx)*DeltaR
C
         Do 30 Ipy = 0, nbpty
            pointb = ymin - rSpher + Dble(Ipy)*DeltaR
C
            Do 40 Ipz = 0, nbptz
               pointc = zmin - rSpher + Dble(Ipz)*DeltaR
C
C
c     IS THIS POINT WITHIN THE VAN DER WAALS ENVELOPE ?
c          Retain1 => false
c     IS THIS POINT BEYOND THE EXCLUSION RADIUS FROM ALL ATOMS ?
c          Retain2 => false
C
               Retain1 = .True.
               Retain2 = .False.
C
               iAtom = 0
50             iAtom = iAtom + 1
               If (iAtom.gt.nAtoms) Goto 60
               If(IsMM(iAtom).ne.0) Goto 50
               RZ = RVDW(IAn(iAtom))
               VDWEnv = RZ*VDWFact
               Distance = (pointa - CO(1,iAtom))**2 +
     &                    (pointb - CO(2,iAtom))**2 +
     &                    (pointc - CO(3,iAtom))**2
               Distance = Sqrt(Distance)
               If(Retain1) Retain1 = Distance .gt. VDWEnv
               Retain2 = Retain2 .or. (Distance .le. (rmax*RZ))
               If(Retain1) Goto 50
C
c     STORE THE POINTS REQUIRED FOR THE FITTING PROCEDURE.
C
60             If (Retain1.and.Retain2) Then
                   nPts = nPts + 1
                   If (Process) Then
                      Grid(1,nPts) = pointa
                      Grid(2,nPts) = pointb
                      Grid(3,nPts) = pointc
                   End If
               End If
40          Continue
30       Continue
20    Continue
C
      Return
      End
