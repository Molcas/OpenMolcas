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
      Subroutine Cavitation(DoDeriv,ToAng,nAt,nS,nTs,GCavP,VMol,TAbs,
     &                      RSolv,Sphere,Tessera,iSphe)
      Implicit Real*8 (A-H,O-Z)
      Logical DoDeriv
C
C     Compute cavitation energy given the cavity.
C
#include "WrkSpc.fh"
      Dimension Sphere(4,*),Tessera(4,*),iSphe(*)
*
      Call GetMem('CavSph'  ,'Allo','Real',ip_CavS,nS      )
      Call GetMem('dCav'    ,'Allo','Real',ip_dCav,3*nAt   )
      Call GetMem('ExpArea' ,'Allo','Real',ip_EA  ,nS      )
      Call GetMem('dExpArea','Allo','Real',ip_dEA ,3*nAt*nS)
      Call FZero(Work(ip_CavS),nS)
      Call FZero(Work(ip_dCav),3*nAt)
      Call FZero(Work(ip_EA),nS)
      Call FZero(Work(ip_dEA),3*nAt*nS)
*
*     Put some quantity in Angstrom
      Call DScal_(nS,ToAng,Sphere(4,1),4)
      Call DScal_(nTs,ToAng*ToAng,Tessera(4,1),4)
*
*     Compute the exposed area for each sphere
      Do iTs = 1, nTs
        iS = iSphe(iTs)
        Work(ip_EA-1+iS) = Work(ip_EA-1+iS) + Tessera(4,iTs)
      End Do
*
      Call Cavitation_(nAt,nS,VMol,TAbs,RSolv,GCavP,Work(ip_CavS),
     &                 Work(ip_dCav),Sphere,Work(ip_EA),Work(ip_dEA),
     &                 DoDeriv)
*
*     Revert to a.u.
      ToAU = 1.0d0 / ToAng
      Call DScal_(nS,ToAU,Sphere(4,1),4)
      Call DScal_(nTs,ToAU*ToAU,Tessera(4,1),4)
*
      Call GetMem('dExpArea','Free','Real',ip_dEA ,3*nAt*nS)
      Call GetMem('ExpArea' ,'Free','Real',ip_EA  ,nS      )
      Call GetMem('dCav'    ,'Free','Real',ip_dCav,3*nAt   )
      Call GetMem('CavSph'  ,'Free','Real',ip_CavS,nS      )
*
      Return
      End
      Subroutine Cavitation_(NAtoms,NSph,VMol,TAbs,RSolv,
     $  GCavP,PCvSph,dCav,Sphere,Ae,dAe,DoDeriv)
      Implicit Real*8 (A-H,O-Z)
      Logical DoDeriv
C
C     Compute cavitation energy given the cavity.
C
      Real*8 Sphere(4,*),Ae(*),dAe(3,NAtoms,*),PCvSph(*),dCav(3,*)
C     Save Zero,One,Two,Three,Four,F4Pt5,F1000,GC,Avgdr,CEKM
      Save Zero,One,Two,Three,Four,F4Pt5,F1000,GC,Avgdr
      Data Zero,One,Two,Three,Four,F4Pt5,F1000/0.0d0,1.0d0,2.0d0,3.0d0,
     $  4.0d0,4.5d0,1.0d3/
C     Data GC,Avgdr,CEKM/1.9865d0,0.60228d0,0.0014389d0/
      Data GC,Avgdr/1.9865d0,0.60228d0/
C
      PI     = Four * ATan(One)
      TPI    = Two  * PI
      FPI    = Two  * TPI
C
C     Cavitation Energy:  R.A. Pierotti, Chem.Rev. 76,717,(1976).
C
C     The cavitation energy is now computed for each sphere of the
C     cavity. The values are scaled for the exposed surface of each
C     sphere and summed irrespective of the number of cavities.
C     The derivative of the cavitation energy wrt nuclear coordinates
C     is also computed. To this end the Pierotti cavitation energy for
C     each sphere is stored in the array PCvSph(ISph).
C
      DENum = Avgdr/VMol
      RT    = GC*TAbs/F1000
      YP    = DENum*FPI*RSolv**3/Three
      YM    = YP/(One-YP)
C
C     Loop on spheres.
C
      GCavP = Zero
      Do 50 ISph = 1, NSph
        RSph = Sphere(4,ISph)
        R    = RSph/RSolv
C
C       Pierotti's cavitation free energy for the single sphere.
C
        G = -log(One-YP) + F4Pt5*YM**2*R**2 + Three*YM*R*(R+One)
        G = RT*G
        PCvSph(ISph) = G
        AExpsd = Ae(ISph)
        ATotal = FPI*RSph*RSph
        Frac   = AExpsd/ATotal
        GCavP  = GCavP + G*Frac
   50 Continue
C
C     Computation of derivatives wrt nuclear coordinates according to
C     the modification of Claverie to the Pierotti approach.
C
C     Loop on atoms, coordinates, tesserae and then spheres.
C
      If(DoDeriv) Then
        Do 70 ISph = 1, NSph
          RSph = Sphere(4,ISph)
          Fact = PCvSph(ISph) / ( FPI * RSph * RSph )
          Do 80 IAtom = 1, NAtoms
            Do 90 Ixyz = 1, 3
              dCav(Ixyz,IAtom) =
     $          dCav(Ixyz,IAtom) + Fact*dAe(Ixyz,IAtom,ISph)
   90       Continue
   80     Continue
   70   Continue
      EndIf
      Return
      End
