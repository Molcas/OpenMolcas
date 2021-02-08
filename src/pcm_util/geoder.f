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
      Subroutine GeoDer(nAt,Cond,nTs,nS,Eps,Sphere,ISphe,NOrd,Tessera,
     $                  Q,DerDM,Grd,DerTes,DerPunt,DerRad,DerCentr)
      Implicit Real*8 (a-h,o-z)
      Dimension Sphere(4,*),ISphe(*),NOrd(*)
      Dimension Tessera(4,*),Q(2,*),DerDM(nTs,*),Grd(3,*)
      Dimension DerTes(nTs,NAt,3),DerPunt(nTs,NAt,3,3)
      Dimension DerRad(nS,NAt,3),DerCentr(nS,NAt,3,3)
      Logical Cond
      Data ToAng /0.52917721067d0/
      Data Zero, Two /0.0d0,2.0d0/
*
*---- Compute the PCM geometric contribution to gradients
*
      Call FZero(Grd,3*nAt)
      Call FZero(DerDM,nTs*nTs)
      Do 1000 IAtom = 1, nAt
        Do 1001 IXYZ = 1, 3
*---- Dielectric model
          If(.not.Cond) then
            Call Over(IAtom,IXYZ,GeoGrd,nAt,nTs,nS,Eps,Sphere,ISphe,
     $                NOrd,Tessera,Q,DerRad,DerCentr)
*---- Conductor model
          ElseIf(Cond) then
            GeoGrd = Zero
            Call DerD(ToAng,IAtom,IXYZ,Tessera,ISphe,DerDM,DerTes,
     $                DerPunt,DerCentr,nTs,nAt,nS)
            Do 2000 iTs = 1, nTs
              Qi = Q(1,iTs) + Q(2,iTs)
              Do 2010 jTs = 1, nTs
                Qj = Q(1,jTs) + Q(2,jTs)
                GeoGrd = GeoGrd + Qi * DerDM(iTs,jTs) * Qj
 2010         Continue
 2000       Continue
          EndIf
          Grd(IXYZ,IAtom) =  GeoGrd / Two
 1001   Continue
 1000 Continue
      Return
      End
*
      Subroutine DERD(ToAng,NSJ,ICOORD,Tessera,ISphe,DerDM,DerTes,
     $                DerPunt,DerCentr,NTs,NAt,NEsf)
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension Tessera(4,*),ISphe(*)
      Dimension DerDM(NTs,*)
      Dimension DerTes(NTs,NAt,3), DerPunt(NTs,NAt,3,3)
      Dimension DerCentr(NEsf,NAt,3,3)
      Data One,Two,Four /1.0d0,2.0d0,4.0d0/
C
      PI  = Four*ATan(One)
      FPI = Four*PI
      ANTOAU = One/ToAng
C
C     Calcola la matrice DERDM, derivata di DMAT (secondo COSMO)
C     rispetto alla coordinata ICOORD della sfera NSJ
C
C     Loop sugli elementi di DMAT
      DO 2000 ITS = 1, NTS
        L = ISPHE(ITS)
        DO 2010 JTS = 1, NTS
          LJ = ISPHE(JTS)
C         Elementi diagonali di DMAT(x)
          IF (ITS.EQ.JTS) THEN
            FAC = - 1.07d0 * Sqrt(FPI) / Two
            DERDM(ITS,ITS) = FAC*DERTES(ITS,NSJ,ICOORD)*ANTOAU /
     $      ( Tessera(4,ITS) * Sqrt(Tessera(4,ITS)) )
          ELSE
C         Elementi fuori diagonale di DMAT(x)
            XIJ = Tessera(1,ITS) - Tessera(1,JTS)
            YIJ = Tessera(2,ITS) - Tessera(2,JTS)
            ZIJ = Tessera(3,ITS) - Tessera(3,JTS)
            DIJ = SQRT( XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ )
            DXIJ=DERPUNT(ITS,NSJ,ICOORD,1)+DERCENTR(L,NSJ,ICOORD,1)
     $          -DERPUNT(JTS,NSJ,ICOORD,1)-DERCENTR(LJ,NSJ,ICOORD,1)
            DYIJ=DERPUNT(ITS,NSJ,ICOORD,2)+DERCENTR(L,NSJ,ICOORD,2)
     $          -DERPUNT(JTS,NSJ,ICOORD,2)-DERCENTR(LJ,NSJ,ICOORD,2)
            DZIJ=DERPUNT(ITS,NSJ,ICOORD,3)+DERCENTR(L,NSJ,ICOORD,3)
     $          -DERPUNT(JTS,NSJ,ICOORD,3)-DERCENTR(LJ,NSJ,ICOORD,3)
            PROD = (XIJ*DXIJ + YIJ*DYIJ + ZIJ*DZIJ) / DIJ**3
            DERDM(ITS,JTS) = - PROD
          ENDIF
 2010 continue
 2000 continue
      RETURN
      END
*
      Subroutine Over(NSJ,ICOORD,GeoGrd,NAt,NTs,NEsf,Eps,Sphere,ISphe,
     $  NOrd,Tessera,Q,DerRad,DerCentr)
      Implicit Real*8(A-H,O-Z)
      Dimension Sphere(4,*),ISphe(*),Tessera(4,*),Q(2,*),NOrd(*)
      Dimension DerRad(NEsf,NAt,3),DerCentr(NEsf,NAt,3,3)
      Data Zero, One, Two, Four /0.0d0,1.0d0,2.0d0,4.0d0/
C
C     Calcola la sovrapposizione delle densita' superficiali sulla
C     porzione di superficie che appartiene alla sfera
C     contenente l'atomo NSJ rispetto al quale stiamo derivando.
C
C     NESFJ e' la sfera che sta attorno all'atomo NSJ: se NSJ non ha
C     nessuna sfera, NESFJ = 0
C
      NESFJ = 0
      DO 2000 I = 1, NESF
        IF( NSJ.EQ.NORD(I) ) NESFJ = I
 2000 continue
*
      SESE=ZERO
      SNSN=ZERO
      SESN=ZERO
      DO 2010 ITS=1,NTS
        L = ISPHE(ITS)
        XNI = - (Sphere(1,L) - Tessera(1,ITS)) / Sphere(4,L)
        YNI = - (Sphere(2,L) - Tessera(2,ITS)) / Sphere(4,L)
        ZNI = - (Sphere(3,L) - Tessera(3,ITS)) / Sphere(4,L)
*
        If(L.EQ.NESFJ) then
         DN = ZERO ! dummy initialize
         If(ICOORD.EQ.1) DN = XNI
         If(ICOORD.EQ.2) DN = YNI
         If(ICOORD.EQ.3) DN = ZNI
        Else
         DCENTN = XNI * DerCentr(L,NSJ,ICOORD,1)
     $          + YNI * DerCentr(L,NSJ,ICOORD,2)
     $          + ZNI * DerCentr(L,NSJ,ICOORD,3)
         DN = DERRAD(L,NSJ,ICOORD) + DCENTN
        EndIf
*
         SNSN = SNSN + DN * Q(1,ITS)**2 / Tessera(4,ITS)
         SESE = SESE + DN * Q(2,ITS)**2 / Tessera(4,ITS)
         SESN = SESN + DN * Q(1,ITS) * Q(2,ITS) / Tessera(4,ITS)
 2010 continue
      PI = Four * ATan(One)
      Fact = Four * PI * Eps / (Eps - One)
      GeoGrd = Fact * (SESE + SNSN + Two*SESN)
      RETURN
      END
