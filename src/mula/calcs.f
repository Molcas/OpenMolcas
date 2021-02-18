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
* Copyright (C) 1994, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine BondStr(R,i1,i2,j,S,NumOfAt,NumInt)
C!
C!  Purpose:
C!    Calculate contribution to the S-matrix due to bond stretching.
C!
C!  Input:
C!    R        : Array of Real*8 real -  contains the
C!               cartesian coordinates of the bond.
C!    i1,i2    : Integer - the number of the atom.
C!    j        : Integer - the number of the internal coordinate.
C!
C!  Output:
C!    S        : Real*8 three dimensional array - the
C!               contributions to S for the parameters specified
C!               in the input.
C!
C!  Uses:
C!    Linalg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 S(3,NumOfAt,NumInt)
      Real*8 R(3)
C!
C!---- Contributions to S.
      SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
      Do k = 1,3
      S(k,i1,j) =-R(k)/SR
      S(k,i2,j) =-S(k,i1,j)
      Call NaNChk(S(k,i1,j),k,i1,j)
      Call NaNChk(S(k,i2,j),k,i1,j)
      End Do
C!
      End
C!
C!-----------------------------------------------------------------------!

      Subroutine NaNChk(X,k,i,j)
      Implicit Real*8 ( a-h,o-z )
      Real*8 X
      Character*16 str1,str2
      Write(str2,'(G16.8)') X
      Call Normalize(str2,str1)
      If(str1(1:3).eq.'NAN') Then
      Write(6,*)' CalcS subroutine produced Not-a-Number!'
      Write(6,*)' Internal coordinate nr j=',j
      Write(6,*)' Atom nr.               i=',i
      Write(6,*)' Component              k=',k
      Write(6,*)' S(k,i,j)=', X
      End If
      Return
      End

C!-----------------------------------------------------------------------!
C!
      Subroutine AngBend(R1,R2,i1,i2,i3,j,S,NumOfAt,NumInt)
C!
C!  Purpose:
C!    Calculate contribution to the S-matrix due to valence
C!    angle bending.
C!
C!  Input:
C!    R1,R2    : Array of Real*8 real - contains the
C!               cartesian coordinates of the bond.
C!    i1,i2,i3 : Integer - the number of the atom.
C!    j        : Integer - the number of the internal coordinate.
C!
C!  Output:
C!    S        : Real*8 three dimensional array - the
C!               contributions to S for the parameters specified
C!               in the input.
C!
C!  Uses:
C!    Linalg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 S(3,NumOfAt,NumInt)
      Real*8 R1(3),R2(3)
      Real*8 NR1(3),NR2(3)
      Real*8 CosTheta,SinTheta,F1,F2,F3
      Real*8 SR1, SR2
C!
C!---- Unit vectors.
      SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
      SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
      Do k = 1,3
      NR1(k) = R1(k)/SR1
      NR2(k) = R2(k)/SR2
      End Do
C!
C!---- Angle Theta.
      Theta = ACos(dDot_(3,NR1,1,NR2,1))
      CosTheta = Cos(Theta)
      SinTheta = Sin(Theta)
C!
C!---- Contributions to S.
      F1 = 1.0d0/(SR1*SinTheta)
      S(1,i1,j) = (CosTheta*NR1(1)-NR2(1))*F1
      S(2,i1,j) = (CosTheta*NR1(2)-NR2(2))*F1
      S(3,i1,j) = (CosTheta*NR1(3)-NR2(3))*F1
C!
      F2 = SR1-SR2*CosTheta
      F3 = SR2-SR1*CosTheta
      F4 = F1/SR2
      S(1,i2,j) = (F2*NR1(1)+F3*NR2(1))*F4
      S(2,i2,j) = (F2*NR1(2)+F3*NR2(2))*F4
      S(3,i2,j) = (F2*NR1(3)+F3*NR2(3))*F4
C!
      F5 = 1.0d0/(SR2*SinTheta)
      S(1,i3,j) = (CosTheta*NR2(1)-NR1(1))*F5
      S(2,i3,j) = (CosTheta*NR2(2)-NR1(2))*F5
      S(3,i3,j) = (CosTheta*NR2(3)-NR1(3))*F5
C!
      Call NaNChk(S(1,i1,j),1,i1,j)
      Call NaNChk(S(2,i1,j),2,i1,j)
      Call NaNChk(S(3,i1,j),3,i1,j)
      Call NaNChk(S(1,i2,j),1,i2,j)
      Call NaNChk(S(2,i2,j),2,i2,j)
      Call NaNChk(S(3,i2,j),3,i2,j)
      Call NaNChk(S(1,i3,j),1,i3,j)
      Call NaNChk(S(2,i3,j),2,i3,j)
      Call NaNChk(S(3,i3,j),3,i3,j)
      End
C!
C!-----------------------------------------------------------------------!
C!-----------------------------------------------------------------------!
C!
      Subroutine LinBend(R1,R2,i1,i2,i3,j,S,NumOfAt,NumInt)
C!
C!  Purpose:
C!    Calculate contribution to the S-matrix due to valence
C!    angle bending for a linear molecule.
C!
C!  Input:
C!    R1,R2    : Array of Real*8 real -  contains the
C!               cartesian coordinates of the bond.
C!    i1,i2,i3 : Integer - the number of the atom.
C!    j        : Integer - the number of the internal coordinate.
C!
C!  Output:
C!    S        : Real*8 three dimensional array - the
C!               contributions to S for the parameters specified
C!               in the input.
C!
C!  Uses:
C!    Linalg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 S(3,NumOfAt,NumInt)
      Real*8 R1(3),R2(3)
C!
C!---- Length of vectors.
      SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
      SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
C!
C!---- Contributions to S.
      F1 = 1.0d0/SR1
      F2 = 1.0d0/SR2
C!
      S(1,i1,  j) = F1
      S(2,i1,  j) = 0.0d0
      S(3,i1,  j) = 0.0d0
C!
      S(1,i2,  j) =-F1-F2
      S(2,i2,  j) = 0.0d0
      S(3,i2,  j) = 0.0d0
C!
      S(1,i3,  j) = F2
      S(2,i3,  j) = 0.0d0
      S(3,i3,  j) = 0.0d0
C!
      S(1,i1,j+1) = 0.0d0
      S(2,i1,j+1) = F1
      S(3,i1,j+1) = 0.0d0
C!
      S(1,i2,j+1) = 0.0d0
      S(2,i2,j+1) =-F1-F2
      S(3,i2,j+1) = 0.0d0
C!
      S(1,i3,j+1) = 0.0d0
      S(2,i3,j+1) = F2
      S(3,i3,j+1) = 0.0d0
C!
      Call NaNChk(S(1,i1,j),1,i1,j)
      Call NaNChk(S(2,i1,j),2,i1,j)
      Call NaNChk(S(3,i1,j),3,i1,j)
      Call NaNChk(S(1,i2,j),1,i2,j)
      Call NaNChk(S(2,i2,j),2,i2,j)
      Call NaNChk(S(3,i2,j),3,i2,j)
      Call NaNChk(S(1,i3,j),1,i3,j)
      Call NaNChk(S(2,i3,j),2,i3,j)
      Call NaNChk(S(3,i3,j),3,i3,j)
      End
C!
C!-----------------------------------------------------------------------!


c!-----------------------------------------------------------------------!
C!
      Subroutine Torsion(R1,R2,R3,i1,i2,i3,i4,j,S,NumOfAt,NumInt)
C!
C!  Purpose:
C!    Calculate contribution to the S-matrix due to torsion.
C!
C!  Input:
C!    R1,R2,R3 : Array of Real*8 real -  contains the
C!               cartesian coordinates of the bond.
C!    i1,i2,
C!    i3,i4    : Integer - the number of the atom.
C!    j        : Integer - the number of the internal coordinate.
C!
C!  Output:
C!    S        : Real*8 three dimensional array - the
C!               contributions to S for the parameters specified
C!               in the input.
C!
C!  Uses:
C!    Linalg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 S(3,NumOfAt,NumInt)
      Real*8 R1(3),R2(3),R3(3)
      Real*8 NR1(3),NR2(3),NR3(3)
C!
C!---- Unit vectors.
      SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
      SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
      SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
      Do k = 1,3
      NR1(k) = R1(k)/SR1
      NR2(k) = R2(k)/SR2
      NR3(k) = R3(k)/SR3
      End Do
C!
C!---- Angles Theta1 and Theta2.
      Theta1 = ACos(-dDot_(3,NR1,1,NR2,1))
      Theta2 = ACos(-dDot_(3,NR2,1,NR3,1))
      CosTheta1 = Cos(Theta1)
      SinTheta1 = Sin(Theta1)
      CosTheta2 = Cos(Theta2)
      SinTheta2 = Sin(Theta2)
C!
C!---- Contributions to S.
      F1 = 1.0d0/(SR1*(SinTheta1)**2)
      F2 = NR1(2)*NR2(3)-NR1(3)*NR2(2)
      F3 = NR1(3)*NR2(1)-NR1(1)*NR2(3)
      F4 = NR1(1)*NR2(2)-NR1(2)*NR2(1)
      S(1,i1,j) =-F2*F1
      S(2,i1,j) =-F3*F1
      S(3,i1,j) =-F4*F1
C!
      F5 = (SR2-SR1*CosTheta1)/(SR1*SR2*(SinTheta1)**2)
      F6 = CosTheta2/(SR2*(SinTheta2)**2)
      F7 = NR3(2)*NR2(3)-NR3(3)*NR2(2)
      F8 = NR3(3)*NR2(1)-NR3(1)*NR2(3)
      F9 = NR3(1)*NR2(2)-NR3(2)*NR2(1)
      S(1,i2,j) = F5*F2+F6*F7
      S(2,i2,j) = F5*F3+F6*F8
      S(3,i2,j) = F5*F4+F6*F9
C!
      F10 = (SR2-SR3*CosTheta2)/(SR2*SR3*SinTheta2**2)
      F11 = CosTheta1/(SR2*(SinTheta1)**2)
      S(1,i3,j) = F10*F7+F11*F2
      S(2,i3,j) = F10*F8+F11*F3
      S(3,i3,j) = F10*F9+F11*F4
C!
      F12 = 1.0d0/(SR3*(SinTheta2)**2)
      S(1,i4,j) =-F7*F12
      S(2,i4,j) =-F8*F12
      S(3,i4,j) =-F9*F12
C!
      Call NaNChk(S(1,i1,j),1,i1,j)
      Call NaNChk(S(2,i1,j),2,i1,j)
      Call NaNChk(S(3,i1,j),3,i1,j)
      Call NaNChk(S(1,i2,j),1,i2,j)
      Call NaNChk(S(2,i2,j),2,i2,j)
      Call NaNChk(S(3,i2,j),3,i2,j)
      Call NaNChk(S(1,i3,j),1,i3,j)
      Call NaNChk(S(2,i3,j),2,i3,j)
      Call NaNChk(S(3,i3,j),3,i3,j)
      Call NaNChk(S(1,i1,j),1,i4,j)
      Call NaNChk(S(2,i1,j),2,i4,j)
      Call NaNChk(S(3,i1,j),3,i4,j)
      End
C!
C!-----------------------------------------------------------------------!
C!-----------------------------------------------------------------------!
C!
      Subroutine OutOfPl(R1,R2,R3,i1,i2,i3,i4,j,S,NumOfAt,NumInt)
C!
C!  Purpose:
C!    Calculate contribution to the S-matrix due to the angle between
C!    a bond and a plane defined by two bonds.
C!
C!  Input:
C!    R1,R2,R3 : Array of Real*8 real -  contains the
C!               cartesian coordinates of the bond.
C!    i1,i2,
C!    i3,i4    : Integer - the number of the atom.
C!    j        : Integer - the number of the internal coordinate.
C!
C!  Output:
C!    S        : Real*8 three dimensional array - the
C!               contributions to S for the parameters specified
C!               in the input.
C!  Uses:
C!    Linalg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 S(3,NumOfAt,NumInt)
      Real*8 R1(3),R2(3),R3(3)
      Real*8 NR1(3),NR2(3),NR3(3),F0(3)
C!
C!---- Unit vectors.
      SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
      SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
      SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
      Do k = 1,3
      NR1(k) = R1(k)/SR1
      NR2(k) = R2(k)/SR2
      NR3(k) = R3(k)/SR3
      End Do
C!
C!---- Calculation of the angle Theta.
      Dot = Ddot_(3,NR2,1,NR3,1)
      Theta = ACos(Dot)
      CosTheta = Cos(Theta)
      SinTheta = Sin(Theta)
C!
C!---- Calculation of the out of plane angle Phi.
      F1 = NR2(2)*NR3(3)-NR2(3)*NR3(2)
      F2 = NR2(3)*NR3(1)-NR2(1)*NR3(3)
      F3 = NR2(1)*NR3(2)-NR2(2)*NR3(1)
      F0(1) = F1
      F0(2) = F2
      F0(3) = F3
      do iv=1,3
      F0(iv) = F0(iv)/SinTheta
      enddo
      Phi = ASin(Ddot_(3,F0,1,NR1,1))
      CosPhi = Cos(Phi)
      TanPhi = Tan(Phi)
C!
C!---- Contributions to S.
      F4 = 1.0d0/(SR1*CosPhi*SinTheta)
      F5 = TanPhi/SR1
      S(1,i1,j) = F1*F4-NR1(1)*F5
      S(2,i1,j) = F2*F4-NR1(2)*F5
      S(3,i1,j) = F3*F4-NR1(3)*F5
C!
      F6 = NR3(2)*NR1(3)-NR3(3)*NR1(2)
      F7 = NR3(3)*NR1(1)-NR3(1)*NR1(3)
      F8 = NR3(1)*NR1(2)-NR3(2)*NR1(1)
      F9 = 1.0d0/(SR2*CosPhi*SinTheta)
      F10 = TanPhi/(SR2*(SinTheta)**2)
      S(1,i2,j) = F6*F9-NR2(1)*F10+CosTheta*F10*NR3(1)
      S(2,i2,j) = F7*F9-NR2(2)*F10+CosTheta*F10*NR3(2)
      S(3,i2,j) = F8*F9-NR2(3)*F10+CosTheta*F10*NR3(3)
C!
      F11 = NR1(2)*NR2(3)-NR1(3)*NR2(2)
      F12 = NR1(3)*NR2(1)-NR1(1)*NR2(3)
      F13 = NR1(1)*NR2(2)-NR1(2)*NR2(1)
      F14 = 1.0d0/(SR3*CosPhi*SinTheta)
      F15 = TanPhi/(SR3*(SinTheta)**2)
      S(1,i3,j) = F11*F14-NR3(1)*F15+CosTheta*F15*NR2(1)
      S(2,i3,j) = F12*F14-NR3(2)*F15+CosTheta*F15*NR2(2)
      S(3,i3,j) = F13*F14-NR3(3)*F15+CosTheta*F15*NR2(3)
C!
      S(1,i4,j) =-S(1,i1,j)-S(1,i2,j)-S(1,i3,j)
      S(2,i4,j) =-S(2,i1,j)-S(2,i2,j)-S(2,i3,j)
      S(3,i4,j) =-S(3,i1,j)-S(3,i2,j)-S(3,i3,j)
C!
      Call NaNChk(S(1,i1,j),1,i1,j)
      Call NaNChk(S(2,i1,j),2,i1,j)
      Call NaNChk(S(3,i1,j),3,i1,j)
      Call NaNChk(S(1,i2,j),1,i2,j)
      Call NaNChk(S(2,i2,j),2,i2,j)
      Call NaNChk(S(3,i2,j),3,i2,j)
      Call NaNChk(S(1,i3,j),1,i3,j)
      Call NaNChk(S(2,i3,j),2,i3,j)
      Call NaNChk(S(3,i3,j),3,i3,j)
      Call NaNChk(S(1,i1,j),1,i4,j)
      Call NaNChk(S(2,i1,j),2,i4,j)
      Call NaNChk(S(3,i1,j),3,i4,j)
      End
C!
C!-----------------------------------------------------------------------!

C!-----------------------------------------------------------------------!
C!
      Subroutine CalcS(AtCoord,InterVec,S,NumInt,NumOfAt)
C!
C!  Purpose:
C!    Calculate the contribution to the S-matrix for each of the
C!    internal coordinates.
C!
C!  Input:
C!    AtCoord  : Two dimensional Real*8 array - contains
C!               the cartesian coordinates of the atoms.
C!    InterVec : Integer array - contains the atoms that are used
C!               in the calculations of each internal coordinate.
C!    S        : Real*8 array filled with zeros.
C!
C!  Output:
C!    S        : Real*8 array - contains the contributions
C!               of each of the internal coordinates.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8  AtCoord (3,NumOfAt)
      Integer InterVec(*)
      Real*8 S(3,NumOfAt,NumInt)
      Real*8 R(3),R1(3),R2(3),R3(3)
C!

C!
C!---- Calculate the contributions to the S-matrix for each internal
C!     coordinate specified in the vector InterVec.
      k = 1
      IntType = InterVec(k)
      nInt = 1
      Do While ( nInt.le.NumInt )
      If(IntType.eq.1) Then
c!---- Bond Stretching.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      do iv=1,3
      R(iv) = (AtCoord(iv,i2)-AtCoord(iv,i1))
      enddo
      Call BondStr(R,i1,i2,nInt,S,NumOfAt,NumInt)
      k = k+3
      EndIf
      If(IntType.eq.2) Then
c!---- Valence Angle Bending.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i3)-AtCoord(iv,i2))
      enddo
      Call AngBend(R1,R2,i1,i2,i3,nInt,S,NumOfAt,NumInt)
      k = k+4
      EndIf
      If(IntType.eq.3) Then
c!---- Linear Valence Angle.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i3)-AtCoord(iv,i2))
      enddo
      Call LinBend(R1,R2,i1,i2,i3,nInt-1,S,NumOfAt,NumInt)
      k = k+4
      EndIf
      If(IntType.eq.4) Then
c!---- Torsion.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      i4 = InterVec(k+4)
      do iv=1,3
      R1(iv) = (AtCoord(iv,i2)-AtCoord(iv,i1))
      R2(iv) = (AtCoord(iv,i3)-AtCoord(iv,i2))
      R3(iv) = (AtCoord(iv,i4)-AtCoord(iv,i3))
      enddo
      Call Torsion(R1,R2,R3,i1,i2,i3,i4,
     &                       nInt,S,NumOfAt,NumInt)
      k = k+5
      EndIf
      If(IntType.eq.5) Then
c!---- Out of Plane Angle Bending.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      i4 = InterVec(k+4)
      do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i4))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i4))
      R3(iv) = (AtCoord(iv,i3)-AtCoord(iv,i4))
      enddo
      Call OutOfPl(R1,R2,R3,i1,i2,i3,i4,nInt,
     &                       S,NumOfAt,NumInt)
      k = k+5
      EndIf
      IntType = InterVec(k)
      nInt = nInt+1
      End Do

C!
      End
C!
C!-----------------------------------------------------------------------!


C!-----------------------------------------------------------------------!
C!
      Subroutine CalcG(G,Mass,S,NumInt,NumOfAt)
C!
C!  Purpose:
C!    Calculate G-matrix.
C!
C!  Input:
C!    Mass     : Real*8 array - the mass of the atoms.
C!    S        : Real*8 three dimensional array - the
C!               contributions from all the internal coordinates.
C!
C!  Output:
C!    G        : Real*8 two dimensional array.
C!
C!  Uses:
C!    Constants
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1994.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 G( NumInt,NumInt )
      Real*8 Mass (NumOfAt)
      Real*8 S(3,NumOfAt,NumInt)
C!
C!
      Do i = 1,NumInt
      Do j = 1,NumInt
      GSum = 0.0d0
      Do k = 1,NumOfAt
      GSum = GSum+(1.0d0/(uToAu*Mass(k)))*
     &                                 (S(1,k,i)*S(1,k,j)+
     &                                  S(2,k,i)*S(2,k,j)+
     &                                  S(3,k,i)*S(3,k,j))
      End Do
      G(i,j) = GSum
      End Do
      End Do
C!
      End
