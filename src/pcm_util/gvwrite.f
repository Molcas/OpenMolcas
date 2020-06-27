************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine GVWrite(Index,NTs,NEsfP,NAt,AtmC,IAt,Coor_Sph,Tessera,
     &                   NVert,Vert,ISphe,Q,ivts,MxVert)

      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension NVert(*),Vert(3,MxVert,*),IVTS(MxVert,*),
     &          Tessera(4,nTs),Q(*)
      Dimension AtmC(3,nAt),IAt(nAt),Coor_Sph(4,*),ISphe(nTs)
C
C     Prepare the input file for GeomView (coloured polyhedra)
C
      Lu = 21
      Lu=IsFreeUnit(Lu)
      If(Index.eq.1) call molcas_open(Lu,'GV.off')
      If(Index.eq.2) call molcas_open(Lu,'GV1.off')
      NumV = 0
      Do 10 i = 1, NTs
  10    NumV = NumV + NVert(i)
      Write(Lu,500)'COFF'
      Write(Lu,1000)NumV, NTs, NumV
      K = 0
      Last = 0
      If(Index.eq.2) then
        qmax=0.0d0
        qmin=0.0d0
        Do 2000 i = 1, NTs
          qt = q(i) / Tessera(4,i)
          if(qt.ge.qmax) qmax = qt
          if(qt.le.qmin) qmin = qt
 2000   continue
        Write(Lu,1100)QMin, QMax
      EndIf
      Do 2010 i = 1, NTs
        If(Index.eq.2) qt = q(i) / Tessera(4,i)
        N=ISphe(i)
        If(N.ne.Last) Write(Lu,1500)N
        Last = N
        If(Index.eq.1) Call Colour(NEsfP,NAt,AtmC,
     +                 IAt,Coor_Sph,N,C1,C2,C3)
        If(Index.eq.2) Call Colchg(i,qt,qmax,qmin,C1,C2,C3)
        Do 2020 j = 1, NVert(i)
          IVTS(j,i) = K
          K = K + 1
          write(Lu,2001)(VERT(jcord,j,i),jcord=1,3),C1,C2,C3,0.75,i
 2020 continue
 2010 continue
      Do 20 i = 1, NTs
  20    write(Lu,3000)NVert(i),(IVTS(j,i),j=1,nvert(i))
      Close(Lu)
 500  Format(1x,a)
1000  Format(3i10)
1100  Format('# Minimum and maximum charge density ',2f12.6)
1500  Format('# Sphere number ',i4)
2001  Format('  ',3f16.9,4f5.2,' # Tess. ',i4)
3000  Format('  ',14i10)
      return
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Subroutine Colour(NesfP,NAt,AtmC,IAt,Coor_Sph,
     +                  N,C1,C2,C3)
      IMPLICIT REAL*8(A-H,O-Z)
      Character*20 Col
      Dimension IAt(*),AtmC(3,*),Coor_Sph(4,*)
C
C     Assign tesserae colours for GeomView:
C     Carbon: green, nitrogen: blue, oxygen: red, hydrogen: light blue,
C     others: violet, added spheres: gray.
C
      Col=' '
      IOut=6
      Delta=1.d-03
      If(N.gt.NESFP)then
       Col = 'Gray'
       Call ColTss(IOut,Col,C1,C2,C3)
       Return
      EndIf
      I = 0
      Do 100 J = 1, NAt
        I = I + 1
        Diff = Sqrt( (AtmC(1,J)-Coor_Sph(1,N))**2
     &             + (AtmC(2,J)-Coor_Sph(2,N))**2
     &             + (AtmC(3,J)-Coor_Sph(3,N))**2 )
        If(Diff.lt.Delta) then
          if(IAt(I).eq.6)then
            Col = 'Green'
          elseif(IAt(I).eq.7)then
            Col = 'Blue'
          elseif(IAt(I).eq.8)then
            Col = 'Red'
          elseif(IAt(I).eq.1)then
            Col = 'Light Blue'
          else
            Col = 'Fuchsia'
          endif
        endif
 100  Continue
      Call ColTss(IOut,Col,C1,C2,C3)
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Subroutine ColTss(IOut,Colour,C1,C2,C3)
      Implicit Real*8(A-H,O-Z)
C
C     Assigne tesserae colours for GeomView:
C
      Character*20 Colour
C
      If(Colour.eq.'White') then
        C1 = 1.0D0
        C2 = 1.0D0
        C3 = 1.0D0
      else if(Colour.eq.'Gray') then
        C1 = 0.66D0
        C2 = 0.66D0
        C3 = 0.66D0
      else if(Colour.eq.'Blue'.or.Colour.eq.'Dark Blue') then
        C1 = 0.0D0
        C2 = 0.0D0
        C3 = 1.0D0
      else if(Colour.eq.'Light Blue') then
        C1 = 0.0D0
        C2 = 1.0D0
        C3 = 1.0D0
      else if(Colour.eq.'Green') then
        C1 = 0.0D0
        C2 = 1.0D0
        C3 = 0.0D0
      else if(Colour.eq.'Yellow') then
        C1 = 1.0D0
        C2 = 1.0D0
        C3 = 0.0D0
      else if(Colour.eq.'Orange') then
        C1 = 1.0D0
        C2 = 0.5D0
        C3 = 0.0D0
      else if(Colour.eq.'Violet') then
        C1 = 0.6D0
        C2 = 0.0D0
        C3 = 1.0D0
      else if(Colour.eq.'Pink'.or.Colour.eq.'Light Red') then
        C1 = 1.0D0
        C2 = 0.5D0
        C3 = 1.0D0
      else if(Colour.eq.'Fuchsia') then
        C1 = 1.0D0
        C2 = 0.0D0
        C3 = 1.0D0
      else if(Colour.eq.'Red'.or.Colour.eq.'Dark Red') then
        C1 = 1.0D0
        C2 = 0.0D0
        C3 = 0.0D0
      else if(Colour.eq.'Black') then
        C1 = 0.0D0
        C2 = 0.0D0
        C3 = 0.0D0
      else
        C1 = 0.0D0  ! dummy assignement
        C2 = 0.0D0  ! dummy assignement
        C3 = 0.0D0  ! dummy assignement
        Write(IOut,'("Unrecognized colour in ColTss")')
        Call Abend()
        endIf
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Subroutine ColChg(i,q,QMAX,QMIN,C1,C2,C3)
      Implicit Real*8(A-H,O-Z)
C
C     Assign tesserae colours for GeomView:
C
      Character*20 Colour
      Save Colour
      Data Colour /'                    '/
C
      Zero = DBLE(0)
      DPos = QMax / 2.d0
      DNeg = QMin / 2.d0
C     Total charge density < 0
      If(q.lt.DNeg) Colour = 'Dark Blue'
      If(q.ge.DNeg.and.q.lt.Zero) Colour = 'Light Blue'
C     Total charge density > 0
      If(q.ge.Zero.and.q.lt.DPos) Colour = 'Pink'
      If(q.ge.DPos) Colour = 'Red'
      Call ColTss(IOut,Colour,C1,C2,C3)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(i)
      End
