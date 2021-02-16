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
* Copyright (C) 1995,1996, Niclas Forsberg                             *
************************************************************************
C!-----------------------------------------------------------------------!
C!
c       Module VibMod
C!
C!  Contains:
C!    VibFreq        (AtCoord,InterVec,Mass,Hess,harmfreq,eigenVec,qMat,
C!                    PED,D3,D4,x_anharm,anharmfreq,max_term)
C!    CalcS          (AtCoord,InterVec,S)
C!    BondStr        (R,i1,i2,j,S)
C!    AngBend        (R1,R2,i1,i2,i3,j,S)
C!    LinBend        (R1,R2,i1,i2,i3,j,S)
C!    Torsion        (R1,R2,R3,i1,i2,i3,i4,j,S)
C!    OutOfPl        (R1,R2,R3,i1,i2,i3,i4,j,S)
C!    CalcG          (G,Mass,S)
C!    Freq           (Hess,G,V,Lambda,B,qMat)
C!    CalcGprime     (Gprime,Mass,xvec,InterVec,AtCoord,NumOfAt)
C!    CalcGdbleprime (Gdbleprime,Mass,xvec,InterVec,AtCoord,NumOfAt)
C!    Anharm         (eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,x)
C!    TransEnergy    (x_anharm,harmfreq,level1,level2)  Result(energy)
C!    AnharmonicFreq (x_anharm,harmfreq,anharmfreq)
C!    Int_to_Cart1   (InterVec,xvec,AtCoord)
C!    Cart_To_Int0   (InterVec,AtCoord,xvec)
C!
C!  Uses:
C!    Constants
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!
C!-----------------------------------------------------------------------!
C!
c         Use IOTools
C!       Implicit Real*8 ( a-h,o-z )
C!
c       Contains

C!-----------------------------------------------------------------------!
C!
      Subroutine VibFreq(AtCoord,xvec,InterVec,Mass,Hess,
     &  G,Gprime,Gdbleprime,
     &    harmfreq,eigenVec,qMat,PED,D3,D4,x_anharm,
     &    anharmfreq,max_term,Cartesian,
     &    nOsc,NumOfAt)
C!
C!  Purpose:
C!    Calculates the vibrational frequencies of a molecule.
C!
C!  Input:
C!    InterVec   : Integer array
C!    Mass       : Real*8 array - masses of the atoms.
C!    xvec       : Real*8 array - geometry of molecule in internal
C!                 coordinates.
C!    Hess       : Real*8 two dimensional array - force constant
C!                 matrix.
C!    D3         : Real*8 three dimensional array - third derivatives
C!                 of potential surface.
C!    D4         : Real*8 four dimensional array - fourth derivatives
C!                 of potential surface.
C!    max_term   : Integer - highest power of term in polynomial fit.
C!    Cartesian  : Logical - If geometry is given in cartesian coordinates,
C!                 then this variable is True.
C!
C!  Output:
C!    AtCoord    : Real*8 two dimensional array - cartesian
C!                 coordinates of the atoms.
C!    harmfreq   : Real*8 array - contains harmonical frequencies.
C!    eigenVec   : Real*8 two dimensional array - contains eigenvectors.
C!    qMat       : Real*8 two dimensional array - cartesian
C!                 displacement vectors.
C!    PED        : Real*8 three dimensional array - potential
C!                 energy distribution.
C!    x_anharm   : Real*8 two dimensional array - anharmonicity
C!                 constants.
C!    anharmfreq : Real*8 array - contains anharmonical frequencies.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
      Integer InterVec(*)
      Real*8 AtCoord (3,NumOfAt)
      Real*8 Mass (NumOfAt)
      Real*8 xvec (nosc)
      Real*8 Hess (nOsc,nOsc)
      Real*8 G (nosc,nosc)
      Real*8 harmfreq (nosc)
      Real*8 anharmfreq (nosc)
      Real*8 eigenVec  (nosc,nosc)
      Real*8 qMat (3*NumOfAt,nOsc)
      Real*8 PED (nosc,nosc,nosc)
      Real*8 D3 (ngdim,ngdim,ngdim)
      Real*8 D4 (ngdim,ngdim,ngdim,ngdim)
      Real*8 x_anharm (nosc,nOsc)
      Real*8 Gprime(ngdim,ngdim,ngdim)
      Real*8 Gdbleprime(ngdim,ngdim,ngdim,ngdim)


      Logical     Cartesian
#include "WrkSpc.fh"
C!
C!---- Initialize.
C!D Write(6,*)' Entered VIBFREQ.'
      NumInt  = nOsc
C!D Write(6,*)' NumInt:',NumInt
C!D Write(6,*)' NumOfAt:',NumOfAt
      Call GetMem('S','Allo','Real',ipS,3*NumOfAt*NumInt)
      Call GetMem('V','Allo','Real',ipV,NumInt*NumInt)
      n=3*NumOfAt
      Call GetMem('B','Allo','Real',ipB,3*NumOfAt*NumInt)

      Call GetMem('Bnew','Allo','Real',ipBnew,3*NumOfAt*NumInt)
      Call GetMem('Lambda','Allo','Real',ipLambda,NumInt)
      call dcopy_(3*NumOfAt*NumInt,0.0d0,0,Work(ipS),1)
C!
C!---- Transform coordinates.
c       xvec=0.0D0
      call dcopy_(nosc,0.0d0,0,xvec,1)
C!D Write(6,*)' VIBFREQ, calling Cart_to_Int0.'
      Call Cart_To_Int0(InterVec,AtCoord,xvec,NumOfAt,NumInt)
C!D Write(6,*)' VIBFREQ, back from Cart_to_Int0.'
C!D Write(6,*)' xvec:'
C!D Write(6,'(5f16.8)') xvec
C!
C!---- Calculate the contributions to the S matrix for each internal
C!     coordinate.
c       S=0.0D0
      call dcopy_(3*NumOfAt*NumInt,0.0d0,0,Work(ipS),1)
C!D Write(6,*)' VIBFREQ, calling CalcS.'
      Call CalcS(AtCoord,InterVec,Work(ipS),NumInt,NumOfAt)
C!D Write(6,*)' VIBFREQ, back from CalcS.'
C!
C!---- Calculate G matrix and first and second derivatives of the G matrix.
C!D Write(6,*)' VIBFREQ, calling CalcG.'
      Call CalcG(G,Mass,Work(ipS),NumInt,NumOfAt)
C!D Write(6,*)' VIBFREQ, back from CalcG.'
c       Gprime = 0.0d0
c       Gdbleprime = 0.0d0
      call dcopy_(ngdim**3,0.0D0,0,GPrime,1)
      call dcopy_(ngdim**4,0.0D0,0,GdblePrime,1)
      If ( max_term.gt.2 ) Then
      dh = 1.0d-3
      Call CalcGprime(Gprime,Mass,xvec,InterVec,AtCoord,
     &    NumOfAt,dh,NumInt)
      dh = 1.0d-2
      Call CalcGdbleprime(Gdbleprime,Mass,xvec,InterVec,
     &    AtCoord,NumOfAt,dh,NumInt)
      End If
C!
C!---- Transform three dimensional array S into two dimensional array B.
      Do j = 1,NumInt
      k = 1
      Do i = 1,NumOfAt
      Work(ipB+k+3*NumOfAt*(j-1)-1) =
     &       Work(ipS+   3*(i-1+NumOfAt*(j-1)))
      Work(ipB+k+3*NumOfAt*(j-1))   =
     &       Work(ipS+1 +3*(i-1+NumOfAt*(j-1)))
      Work(ipB+k+3*NumOfAt*(j-1)+1) =
     &       Work(ipS+2 +3*(i-1+NumOfAt*(j-1)))
      k = k+3
      End Do
      End Do
C!---- Transform three dimensional array S into two dimensional array Bnew.
      Do j = 1,NumInt
      k = 1
      Do i = 1,NumOfAt
      Work(ipBnew+k  +3*NumOfAt*(j-1)-1) =
     &       Work(ipS+   3*(i-1+NumOfAt*(j-1)))/Sqrt(uToAU*Mass(i))
      Work(ipBnew+k  +3*NumOfAt*(j-1))   =
     &       Work(ipS+1 +3*(i-1+NumOfAt*(j-1)))/Sqrt(uToAU*Mass(i))
      Work(ipBnew+k  +3*NumOfAt*(j-1)+1) =
     &       Work(ipS+2 +3*(i-1+NumOfAt*(j-1)))/Sqrt(uToAU*Mass(i))
      k = k+3
      End Do
      End Do
C!
C!---- Given Hess and G, calculate the eigenvalues and eigenvectors of G*Hess.
C!D Write(6,*)' VIBFREQ, calling Freq.'
      Call Freq(Hess,G,Work(ipV),Work(ipLambda),Work(ipB),
     &  Work(ipBnew),qMat,nOsc,NumOfAt)
C!D Write(6,*)' VIBFREQ, back from Freq.'
C!D Write(6,*)' Lambda:'
C!D Write(6,'(5f16.8)') Lambda
C!
C!---- Calculate harmonic frequencies.

      do iv=1,nOsc
      harmfreq(iv) = sqrt(abs(Work(ipLambda+iv-1)))
      enddo
C!D Write(6,*)' harmfreq:'
C!D Write(6,'(5f16.8)') harmfreq
c       eigenVec = V
      call dcopy_(nOsc*nOsc,Work(ipV),1,eigenVec,1)
C!
C!---- Anharmonicity calculations (if we have third and possibly fourth
C!     derivatives). First calculation of the anharmonicity constants and
C!     then calculation of the fundamental frequencies.
c       x_anharm = 0.0d0
      call dcopy_(nosc**2,0.0d0,0,x_anharm,1)
      If ( max_term.gt.2 ) Then
      Call GetMem('C1','Allo','Real',ipC1,nOsc*nOsc)
      Call GetMem('Temp','Allo','Real',ipTemp,nOsc*nOsc)
      Call GetMem('V3','Allo','Real',ipV3,nOsc*nOsc*nOsc)
      Call GetMem('T3','Allo','Real',ipT3,nOsc*nOsc*nOsc)
      Call GetMem('V4','Allo','Real',ipV4,nOsc*nOsc*nOsc*nOsc)
      Call GetMem('T4','Allo','Real',ipT4,nOsc*nOsc*nOsc*nOsc)
      Call Anharm(eigenVec,harmfreq,D3,D4,Gprime,Gdbleprime,
     &    x_anharm,max_term,nOsc,
     &    Work(ipC1),Work(ipTemp),Work(ipV3),Work(ipT3),
     &  Work(ipV4),Work(ipT4))
      Call GetMem('C1','Free','Real',ipC1,nOsc*nOsc)
      Call GetMem('Temp','Free','Real',ipTemp,nOsc*nOsc)
      Call GetMem('V3','Free','Real',ipV3,nOsc*nOsc*nOsc)
      Call GetMem('T3','Free','Real',ipT3,nOsc*nOsc*nOsc)
      Call GetMem('V4','Free','Real',ipV4,nOsc*nOsc*nOsc*nOsc)
      Call GetMem('T4','Free','Real',ipT4,nOsc*nOsc*nOsc*nOsc)

      Call AnharmonicFreq(x_anharm,harmfreq,anharmfreq,nOsc)
      End If
C!
C!---- Calculate potential energy distribution.
      Call PotDist(Hess,Work(ipV),Work(ipLambda),PED,NumInt,nOsc)

C!---- Free memory space of S, B, Bnew, G and V.
      Call GetMem('S','Free','Real',ipS,3*NumOfAt*NumInt)
      Call GetMem('B','Free','Real',ipB,3*NumOfAt*NumInt)
      Call GetMem('V','Free','Real',ipV,NumInt*NumInt)
      Call GetMem('Lambda','Free','Real',ipLambda,NumInt)
      Call GetMem('Bnew','Free','Real',ipBnew,3*NumOfAt*NumInt)
C!
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Cartesian)
      End
