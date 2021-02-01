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
      Subroutine SupSym(FrcCrt,nsAtom,Coor,nSupSy,Idntcl,iAtom)
      use Slapaf_Info, only: dMass, Smmtrc
      Implicit Real*8 (a-h,o-z)
************************************************************************
*                                                                      *
*     Object: To constrain forces to higher symmetry.                  *
*                                                                      *
************************************************************************
#include "real.fh"
*
*           --- global arrays ---
*
      Real*8 FrcCrt(3,nsAtom), Coor(3,nsAtom)
      Integer   Idntcl(nSupSy), iAtom(nsAtom)
*
*         --- local arrays ---
*
      Real*8 cMass(3)
      Real*8 E(3)

      Call CofMss(Coor,nsAtom,cMass)
*
*
*     Loop over groups of centers which are identical.
*
      K = 0
      Do 10 I = 1, nSupSy
         FORCE = Zero
*
*        loop over centers which are identical.
*
         IDIV = 0
         Do 15 J = 1, Idntcl(I)
            K = K + 1
*           Get the number of this center
            ICEN = iAtom(K)
*
*           Get an unit vector in the direction from center of mass to
*           the ICEN'th center.
*
            E2 = Zero
            Do 20 L = 1, 3
               E(L) = Coor(L,ICEN) - cMass(L)
               E2 = E2 + E(L)*E(L)
20          Continue
            ABSE = SQRT(E2)
*           Normalize the vector.
            Do 25 L = 1, 3
               E(L) = E(L) / ABSE
25          Continue
*
*           Project the force of the ICEN'th center on this unit vector
*
            PROJ = DDot_(3,E,1,FrcCrt(1,ICEN),1)
*
*           Sum forces which should be equal
            FORCE = FORCE + PROJ * DBLE(iDeg(Coor(1,iCen)))
            IDIV = IDIV + iDeg(Coor(1,iCen))
*
*           Save the unit vector of the new force of this center
*
            Do 30 L = 1, 3
               FrcCrt(L,ICEN) = E(L)
30          Continue
15       Continue
*
*        Get the new force
         FORCE = FORCE / DBLE(IDIV)
         K = K - Idntcl(I)
*
*        Symmetrize the forces
*
         Do 35 J = 1, Idntcl(I)
            K = K + 1
            ICEN = iAtom(K)
*
            Do 40 L = 1, 3
               FrcCrt(L,ICEN) = FrcCrt(L,ICEN) * FORCE
40          Continue
35       Continue
10    Continue
*
      Contains
      Subroutine CofMss(Coor,nsAtom,cMass)
************************************************************************
*     Object: To calculate the molecular mass, the center of mass and  *
*             move the coordinates so origo is the center of mass.     *
************************************************************************
      Real*8 COOR(3,nsAtom), cMass(3)
      Integer i, j
*
*     Calculate the molecular mass.
*
      TMass = Zero
      Do I = 1, nsAtom
         TMass = TMass + dMass(I) * DBLE(iDeg(Coor(1,i)))
      End Do
      iCOM=-1
      If (TMass.ge.1.D99) Then
         Do i = 1, nsAtom
            If (dMass(i).eq.1.D99) Then
               iCOM=i
               Go To 99
            End If
         End Do
      End If
 99   Continue
*
*     calculate the center of mass
*
      cMass(:)=Zero
*-----Loop over the unique centers
      Do i = 1, nsAtom
         Do j = 1, 3
*-----------Add contribution
            If (Smmtrc(j,i)) cMass(j) = cMass(j) +
     &         dMass(i) *  Coor(j,i) *
     &         DBLE(iDeg(Coor(1,i)))
         End Do
      End Do
*
      Do i = 1, 3
         cMass(i) = cMass(i) / TMass
      End Do
      If (iCOM.ge.1.and.iCOM.le.nsAtom) cMass(:)=Coor(:,iCom)
*
#ifdef _DEBUGPRINT_
      If (LWrite) Write(6,100) (cMass(i),i=1,3), TMass
 100  FORMAT(//,' Center of Mass (Bohr) ',3F10.5,/,
     &          ' Molecular Mass   (au) ',1F15.5)
#endif
#ifdef _DO_NOT_
*
*     translate the center of mass to origo
*
      Do i = 1, nsAtom
         Do j = 1, 3
            Coor(j,i) = Coor(j,i) - cMass(j)
         End Do
      End Doa
#endif
*
      Return
      End Subroutine CofMss

      End Subroutine SupSym
