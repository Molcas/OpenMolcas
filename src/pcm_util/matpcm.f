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
      Subroutine MatPCM(NTs,Eps,Conductor,ISphe,Coor_Sph,Tessera,
     &                  DMat,SMat,SDMat,TMat,RMat)
      Implicit Real*8 (A-H,O-Z)
      Logical Conductor
*  Compute PCM matrix with the formalism in
*  M. Cossi, N. Rega, G. Scalmani, V. Barone JCP in press;
*  D.M. Chipman JCP 112, 5558 (2000).
*  Solvation charges are defined through
*  Tq = RV
*  where V is the solute electrostatic potential. Here T^-1*R is computed
*  and finally returned in DMat.
      Dimension ISphe(*),Tessera(4,nTs)
      Dimension Coor_Sph(4,*)
      Dimension SMat(NTs,*),TMat(NTs,*),RMat(NTs,*),DMat(NTs,*)
      Dimension SDMat(NTs,*)
      Data Zero, One, Two, Four /0.0d0, 1.0d0, 2.0d0, 4.0d0/
      Data PotFac /1.0694d0/
      PI  = Four*ATan(One)
      TPI = Two * PI
      FPI = Two * TPI
      If(Conductor) goto 100
*
* Dielectric model:
*
* S and D* matrices
      call dcopy_(nTs*nTs,[Zero],0,DMat,1)
      Do 1000 ITs = 1, NTs
        XI = Tessera(1,iTs)
        YI = Tessera(2,iTs)
        ZI = Tessera(3,iTs)
        LI = ISphe(ITs)
        XNI = (XI - Coor_Sph(1,LI)) / Coor_Sph(4,LI)
        YNI = (YI - Coor_Sph(2,LI)) / Coor_Sph(4,LI)
        ZNI = (ZI - Coor_Sph(3,LI)) / Coor_Sph(4,LI)
        SMat(ITs,ITs) = PotFac * Sqrt(FPI / Tessera(4,ITs))
        DMat(ITs,ITs) = DMat(ITs,ITs) - TPI / Tessera(4,ITs)
        Do 1001 JTs = 1, NTs
          If(JTs.eq.ITs) goto 1001
          XJ = Tessera(1,jTs)
          YJ = Tessera(2,jTs)
          ZJ = Tessera(3,jTs)
          RIJ = Sqrt( (XI - XJ)**2 + (YI - YJ)**2 + (ZI - ZJ)**2 )
          SMat(ITs,JTs) = One / RIJ
          Prod = (XI-XJ) * XNI + (YI-YJ) * YNI + (ZI-ZJ) * ZNI
          DMat(ITs,JTs) = - Prod / RIJ**3
          DMat(JTs,JTs) = DMat(JTs,JTs)
     &                  - DMat(ITs,JTs)*Tessera(4,ITs) / Tessera(4,JTs)
 1001   Continue
 1000 Continue
*
* S*A*D matrix
      call dcopy_(nTs*nTs,[Zero],0,SDMat,1)
      Do 1500 ITs = 1, NTs
        Do 1501 JTs = 1, NTs
          Do 1502 KTs = 1, NTs
            SDMat(ITs,JTs) = SDMat(ITs,JTs) +
     &      SMat(ITs,KTs) * Tessera(4,KTs) * DMat(KTs,JTs)
 1502     Continue
 1501   Continue
 1500 Continue
*
* The charges are defined as
* q = T-1 R V,         T = f(e)*S - SAD / 2p
*
* T and R matrices
      Fac = (Eps + One) / (Eps - One)
      Do 2000 ITs = 1, NTs
        Rad = Coor_Sph(4,ISphe(ITs))
        TMat(ITs,ITs) = Fac * SMat(ITs,ITs) - SDMat(ITs,ITs) / TPI
        RMat(ITs,ITs) = - One + DMat(ITs,ITs) * Tessera(4,ITs) / TPI
        Do 2001 JTs = 1, NTs
          If(JTs.eq.ITs) goto 2001
          TMat(ITs,JTs) = Fac * SMat(ITs,JTs) - SDMat(ITs,JTs) / TPI
          RMat(ITs,JTs) = Tessera(4,JTs) * DMat(JTs,ITs) / TPI
 2001   Continue
 2000 Continue
*
* Invert T matrix
*
      If (Eps.gt.One) Then
         Call MatInvert(TMat,nTs)
      Else
         Call FZero(TMat,nTs**2)
      End If
*
* Form T^-1 * R and store it in D
*
      Call DGEMM_('N','N',
     &            nTs,nTs,nTs,
     &            1.0d0,TMat,nTs,
     &            RMat,nTs,
     &            0.0d0,DMat,nTs)
      Return
  100 Continue
*
* Conductor model
*
* S matrix
      EpsFac = Eps / (Eps - One)
      call dcopy_(nTs*nTs,[Zero],0,SMat,1)
      Do 1010 ITs = 1, NTs
        XI = Tessera(1,iTs)
        YI = Tessera(2,iTs)
        ZI = Tessera(3,iTs)
        SMat(ITs,ITs) = - PotFac * EpsFac * Sqrt(FPI / Tessera(4,ITs))
        Do 1011 JTs = 1, ITs-1
          XJ = Tessera(1,jTs)
          YJ = Tessera(2,jTs)
          ZJ = Tessera(3,jTs)
          RIJ = Sqrt( (XI - XJ)**2 + (YI - YJ)**2 + (ZI - ZJ)**2 )
          SMat(ITs,JTs) = - EpsFac * One / RIJ
          SMat(JTs,ITs) = SMat(ITs,JTs)
 1011   Continue
 1010 Continue
*
* Invert S matrix and store it in D
*
      If (Eps.gt.One) Then
         Call MatInvert(SMat,nTs)
         call dcopy_(nTs*nTs,SMat,1,DMat,1)
      Else
         Call FZero(DMat,nTs**2)
      End If
      Return
      End
