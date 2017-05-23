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
c     Compute F3 using cumulant reconstruction except for G3-dependent terms
      Function CU4F3H (NAC,E,ES,G1,G2,F1,F2,iP,iQ,jP,jQ,kP,kQ)
      Implicit Real*8 (A-H,O-Z)
      Dimension E(NAC)
      Dimension G1(NAC,NAC), G2(NAC,NAC,NAC,NAC)
      Dimension F1(NAC,NAC), F2(NAC,NAC,NAC,NAC)
c
      CU4F3H=0.0D0
c
c Compute these 3PDM contributions in other place
c
c     // G(I,i)G(JKL,jkl)[4]
c     CU4F3H=CU4F3H
c    & + G1(lT,lT)*G3(iP,iQ,jP,jQ,kP,kQ)
c     // G(I,j)G(JKL,ikl)[12]
c     CU4F3H=CU4F3H
c    & - 0.5D0*G1(iP,lT)*G3(lT,iQ,jP,jQ,kP,kQ)
c    & - 0.5D0*G1(jP,lT)*G3(lT,jQ,iP,iQ,kP,kQ)
c    & - 0.5D0*G1(kP,lT)*G3(lT,kQ,iP,iQ,jP,jQ)
c    & - 0.5D0*G1(lT,iQ)*G3(iP,lT,jP,jQ,kP,kQ)
c    & - 0.5D0*G1(lT,jQ)*G3(jP,lT,iP,iQ,kP,kQ)
c    & - 0.5D0*G1(lT,kQ)*G3(kP,lT,iP,iQ,jP,jQ)
c

c     // G(I,i)G(JKL,jkl)[4]
      CU4F3H=CU4F3H
     & + G1(iP,iQ)*F2(jP,jQ,kP,kQ)
     & + G1(jP,jQ)*F2(iP,iQ,kP,kQ)
     & + G1(kP,kQ)*F2(iP,iQ,jP,jQ)
c     // G(I,j)G(JKL,ikl)[12]
      CU4F3H=CU4F3H
     & - 0.5D0*G1(iP,jQ)*F2(jP,iQ,kP,kQ)
     & - 0.5D0*G1(iP,kQ)*F2(kP,iQ,jP,jQ)
     & - 0.5D0*G1(jP,iQ)*F2(iP,jQ,kP,kQ)
     & - 0.5D0*G1(jP,kQ)*F2(kP,jQ,iP,iQ)
     & - 0.5D0*G1(kP,iQ)*F2(iP,kQ,jP,jQ)
     & - 0.5D0*G1(kP,jQ)*F2(jP,kQ,iP,iQ)
c     // G(IJ,ij)G(KL,kl)[3]
      CU4F3H=CU4F3H
     & + G2(iP,iQ,jP,jQ)*F1(kP,kQ)
     & + G2(iP,iQ,kP,kQ)*F1(jP,jQ)
     & + F1(iP,iQ)*G2(jP,jQ,kP,kQ)
c     // G(IJ,ik)G(KL,jl)[12]
      CU4F3H=CU4F3H
     & - 0.5D0*G2(iP,iQ,jP,kQ)*F1(kP,jQ)
     & - 0.5D0*G2(jP,jQ,iP,kQ)*F1(kP,iQ)
     & - 0.5D0*G2(kP,kQ,iP,jQ)*F1(jP,iQ)
     & - 0.5D0*F1(iP,jQ)*G2(jP,iQ,kP,kQ)
     & - 0.5D0*F1(iP,kQ)*G2(kP,iQ,jP,jQ)
     & - 0.5D0*F1(jP,kQ)*G2(kP,jQ,iP,iQ)
c     // G(I,i)G(J,j)G(KL,kl)[6]
      CU4F3H=CU4F3H
     & - 2.0D0*G1(iP,iQ)*G1(jP,jQ)*F1(kP,kQ)
     & - 2.0D0*G1(iP,iQ)*G1(kP,kQ)*F1(jP,jQ)
     & - 2.0D0*G1(jP,jQ)*G1(kP,kQ)*F1(iP,iQ)
     & - 2.0D0*G1(iP,iQ)*ES*G2(jP,jQ,kP,kQ)
     & - 2.0D0*G1(jP,jQ)*ES*G2(iP,iQ,kP,kQ)
     & - 2.0D0*G1(kP,kQ)*ES*G2(iP,iQ,jP,jQ)
c     // G(I,j)G(J,i)G(KL,kl)[6]
      CU4F3H=CU4F3H
     & + G1(iP,jQ)*G1(jP,iQ)*F1(kP,kQ)
     & + G1(iP,kQ)*G1(kP,iQ)*F1(jP,jQ)
     & + G1(jP,kQ)*G1(kP,jQ)*F1(iP,iQ)
c     // G(I,i)G(J,k)G(KL,jl)[24]
      CU4F3H=CU4F3H
     & + G1(iP,iQ)*G1(jP,kQ)*F1(kP,jQ)
     & + G1(iP,iQ)*G1(kP,jQ)*F1(jP,kQ)
     & + G1(jP,jQ)*G1(iP,kQ)*F1(kP,iQ)
     & + G1(jP,jQ)*G1(kP,iQ)*F1(iP,kQ)
     & + G1(kP,kQ)*G1(iP,jQ)*F1(jP,iQ)
     & + G1(kP,kQ)*G1(jP,iQ)*F1(iP,jQ)
     & + ES*G1(iP,jQ)*G2(jP,iQ,kP,kQ)
     & + ES*G1(jP,iQ)*G2(iP,jQ,kP,kQ)
     & + ES*G1(iP,kQ)*G2(kP,iQ,jP,jQ)
     & + ES*G1(kP,iQ)*G2(iP,kQ,jP,jQ)
     & + ES*G1(jP,kQ)*G2(kP,jQ,iP,iQ)
     & + ES*G1(kP,jQ)*G2(jP,kQ,iP,iQ)
c     // G(I,j)G(J,k)G(KL,il)[24]
      CU4F3H=CU4F3H
     & - 0.5D0*G1(iP,jQ)*G1(jP,kQ)*F1(kP,iQ)
     & - 0.5D0*G1(iP,kQ)*G1(kP,jQ)*F1(jP,iQ)
     & - 0.5D0*G1(jP,iQ)*G1(iP,kQ)*F1(kP,jQ)
     & - 0.5D0*G1(jP,kQ)*G1(kP,iQ)*F1(iP,jQ)
     & - 0.5D0*G1(kP,iQ)*G1(iP,jQ)*F1(jP,kQ)
     & - 0.5D0*G1(kP,jQ)*G1(jP,iQ)*F1(iP,kQ)
c     // G(I,i)G(J,j)G(K,k)G(L,l)[1]
      CU4F3H=CU4F3H
     & + 6.0D0*G1(iP,iQ)*G1(jP,jQ)*G1(kP,kQ)*ES
c     // G(I,j)G(J,i)G(K,k)G(L,l)[6]
      CU4F3H=CU4F3H
     & - 3.0D0*G1(iP,jQ)*G1(jP,iQ)*G1(kP,kQ)*ES
     & - 3.0D0*G1(iP,kQ)*G1(kP,iQ)*G1(jP,jQ)*ES
     & - 3.0D0*G1(jP,kQ)*G1(kP,jQ)*G1(iP,iQ)*ES
c     // G(I,j)G(J,k)G(K,i)G(L,l)[8]
      CU4F3H=CU4F3H
     & + 1.5D0*G1(iP,jQ)*G1(jP,kQ)*G1(kP,iQ)*ES
     & + 1.5D0*G1(iP,kQ)*G1(kP,jQ)*G1(jP,iQ)*ES

c Loop over lT
      DO 900 lT=1,NAC
c     // G(IJ,ik)G(KL,jl)[12]
      CU4F3H=CU4F3H
     & - 0.5D0*G2(iP,iQ,jP,lT)*G2(lT,jQ,kP,kQ)*E(lT)
     & - 0.5D0*G2(iP,iQ,kP,lT)*G2(lT,kQ,jP,jQ)*E(lT)
     & - 0.5D0*G2(jP,jQ,iP,lT)*G2(lT,iQ,kP,kQ)*E(lT)
     & - 0.5D0*G2(jP,jQ,kP,lT)*G2(lT,kQ,iP,iQ)*E(lT)
     & - 0.5D0*G2(kP,kQ,iP,lT)*G2(lT,iQ,jP,jQ)*E(lT)
     & - 0.5D0*G2(kP,kQ,jP,lT)*G2(lT,jQ,iP,iQ)*E(lT)
c     // G(IJ,kl)G(KL,ij)[3], G(IJ,lk)G(KL,ji)[3]
c     // G(IJ,kl)G(KL,ji)[3], G(IJ,lk)G(KL,ij)[3]
      CU4F3H=CU4F3H
     & + (
     &   + 2.0D0*G2(iP,kQ,jP,lT)*G2(kP,iQ,lT,jQ)*E(lT)
     &   + 2.0D0*G2(iP,jQ,kP,lT)*G2(jP,iQ,lT,kQ)*E(lT)
     &   + 2.0D0*G2(iP,jQ,lT,kQ)*G2(jP,iQ,kP,lT)*E(lT)
     &   + 2.0D0*G2(iP,lT,jP,kQ)*G2(kP,jQ,lT,iQ)*E(lT)
     &   + 2.0D0*G2(iP,lT,kP,jQ)*G2(jP,kQ,lT,iQ)*E(lT)
     &   + 2.0D0*G2(iP,kQ,lT,jQ)*G2(jP,lT,kP,iQ)*E(lT)
     &   + G2(iP,kQ,jP,lT)*G2(kP,jQ,lT,iQ)*E(lT)
     &   + G2(iP,jQ,kP,lT)*G2(jP,kQ,lT,iQ)*E(lT)
     &   + G2(iP,jQ,lT,kQ)*G2(jP,lT,kP,iQ)*E(lT)
     &   + G2(iP,lT,jP,kQ)*G2(kP,iQ,lT,jQ)*E(lT)
     &   + G2(iP,lT,kP,jQ)*G2(jP,iQ,lT,kQ)*E(lT)
     &   + G2(iP,kQ,lT,jQ)*G2(jP,iQ,kP,lT)*E(lT)
     & )/6.0D0
c     // G(I,j)G(J,i)G(KL,kl)[6]
      CU4F3H=CU4F3H
     & + G1(iP,lT)*G1(lT,iQ)*G2(jP,jQ,kP,kQ)*E(lT)
     & + G1(jP,lT)*G1(lT,jQ)*G2(iP,iQ,kP,kQ)*E(lT)
     & + G1(kP,lT)*G1(lT,kQ)*G2(iP,iQ,jP,jQ)*E(lT)
c     // G(I,i)G(J,k)G(KL,jl)[24]
      CU4F3H=CU4F3H
     & + G1(iP,iQ)*G1(jP,lT)*G2(lT,jQ,kP,kQ)*E(lT)
     & + G1(iP,iQ)*G1(kP,lT)*G2(lT,kQ,jP,jQ)*E(lT)
     & + G1(jP,jQ)*G1(iP,lT)*G2(lT,iQ,kP,kQ)*E(lT)
     & + G1(jP,jQ)*G1(kP,lT)*G2(lT,kQ,iP,iQ)*E(lT)
     & + G1(kP,kQ)*G1(iP,lT)*G2(lT,iQ,jP,jQ)*E(lT)
     & + G1(kP,kQ)*G1(jP,lT)*G2(lT,jQ,iP,iQ)*E(lT)
     & + G1(iP,iQ)*G1(lT,jQ)*G2(jP,lT,kP,kQ)*E(lT)
     & + G1(iP,iQ)*G1(lT,kQ)*G2(kP,lT,jP,jQ)*E(lT)
     & + G1(jP,jQ)*G1(lT,iQ)*G2(iP,lT,kP,kQ)*E(lT)
     & + G1(jP,jQ)*G1(lT,kQ)*G2(kP,lT,iP,iQ)*E(lT)
     & + G1(kP,kQ)*G1(lT,iQ)*G2(iP,lT,jP,jQ)*E(lT)
     & + G1(kP,kQ)*G1(lT,jQ)*G2(jP,lT,iP,iQ)*E(lT)
c     // G(I,j)G(J,k)G(KL,il)[24]
      CU4F3H=CU4F3H
     & - 0.5D0*G1(iP,jQ)*G1(jP,lT)*G2(lT,iQ,kP,kQ)*E(lT)
     & - 0.5D0*G1(iP,kQ)*G1(kP,lT)*G2(lT,iQ,jP,jQ)*E(lT)
     & - 0.5D0*G1(jP,iQ)*G1(iP,lT)*G2(lT,jQ,kP,kQ)*E(lT)
     & - 0.5D0*G1(jP,kQ)*G1(kP,lT)*G2(lT,jQ,iP,iQ)*E(lT)
     & - 0.5D0*G1(kP,iQ)*G1(iP,lT)*G2(lT,kQ,jP,jQ)*E(lT)
     & - 0.5D0*G1(kP,jQ)*G1(jP,lT)*G2(lT,kQ,iP,iQ)*E(lT)
     & - 0.5D0*G1(iP,lT)*G1(lT,jQ)*G2(jP,iQ,kP,kQ)*E(lT)
     & - 0.5D0*G1(iP,lT)*G1(lT,kQ)*G2(kP,iQ,jP,jQ)*E(lT)
     & - 0.5D0*G1(jP,lT)*G1(lT,iQ)*G2(iP,jQ,kP,kQ)*E(lT)
     & - 0.5D0*G1(jP,lT)*G1(lT,kQ)*G2(kP,jQ,iP,iQ)*E(lT)
     & - 0.5D0*G1(kP,lT)*G1(lT,iQ)*G2(iP,kQ,jP,jQ)*E(lT)
     & - 0.5D0*G1(kP,lT)*G1(lT,jQ)*G2(jP,kQ,iP,iQ)*E(lT)
     & - 0.5D0*G1(lT,iQ)*G1(iP,jQ)*G2(jP,lT,kP,kQ)*E(lT)
     & - 0.5D0*G1(lT,iQ)*G1(iP,kQ)*G2(kP,lT,jP,jQ)*E(lT)
     & - 0.5D0*G1(lT,jQ)*G1(jP,iQ)*G2(iP,lT,kP,kQ)*E(lT)
     & - 0.5D0*G1(lT,jQ)*G1(jP,kQ)*G2(kP,lT,iP,iQ)*E(lT)
     & - 0.5D0*G1(lT,kQ)*G1(kP,iQ)*G2(iP,lT,jP,jQ)*E(lT)
     & - 0.5D0*G1(lT,kQ)*G1(kP,jQ)*G2(jP,lT,iP,iQ)*E(lT)
c     // G(I,k)G(J,l)G(KL,ij)[12]
      CU4F3H=CU4F3H
     & - 0.5D0*G1(iP,kQ)*G1(jP,lT)*G2(kP,iQ,lT,jQ)*E(lT)
     & - 0.5D0*G1(iP,jQ)*G1(kP,lT)*G2(jP,iQ,lT,kQ)*E(lT)
     & - 0.5D0*G1(jP,iQ)*G1(kP,lT)*G2(iP,jQ,lT,kQ)*E(lT)
     & - 0.5D0*G1(iP,lT)*G1(jP,kQ)*G2(lT,iQ,kP,jQ)*E(lT)
     & - 0.5D0*G1(iP,lT)*G1(kP,jQ)*G2(lT,iQ,jP,kQ)*E(lT)
     & - 0.5D0*G1(jP,lT)*G1(kP,iQ)*G2(lT,jQ,iP,kQ)*E(lT)
     & - 0.5D0*G1(iP,jQ)*G1(lT,kQ)*G2(jP,iQ,kP,lT)*E(lT)
     & - 0.5D0*G1(iP,kQ)*G1(lT,jQ)*G2(kP,iQ,jP,lT)*E(lT)
     & - 0.5D0*G1(jP,iQ)*G1(lT,kQ)*G2(iP,jQ,kP,lT)*E(lT)
     & - 0.5D0*G1(jP,kQ)*G1(lT,iQ)*G2(kP,jQ,iP,lT)*E(lT)
     & - 0.5D0*G1(kP,iQ)*G1(lT,jQ)*G2(iP,kQ,jP,lT)*E(lT)
     & - 0.5D0*G1(kP,jQ)*G1(lT,iQ)*G2(jP,kQ,iP,lT)*E(lT)
c     // G(I,j)G(J,i)G(K,k)G(L,l)[6]
      CU4F3H=CU4F3H
     & - 3.0D0*G1(iP,lT)*G1(lT,iQ)*G1(jP,jQ)*G1(kP,kQ)*E(lT)
     & - 3.0D0*G1(jP,lT)*G1(lT,jQ)*G1(iP,iQ)*G1(kP,kQ)*E(lT)
     & - 3.0D0*G1(kP,lT)*G1(lT,kQ)*G1(iP,iQ)*G1(jP,jQ)*E(lT)
c     // G(I,j)G(J,k)G(K,i)G(L,l)[8]
      CU4F3H=CU4F3H
     & + 1.5D0*G1(iP,jQ)*G1(jP,lT)*G1(lT,iQ)*G1(kP,kQ)*E(lT)
     & + 1.5D0*G1(iP,kQ)*G1(kP,lT)*G1(lT,iQ)*G1(jP,jQ)*E(lT)
     & + 1.5D0*G1(jP,kQ)*G1(kP,lT)*G1(lT,jQ)*G1(iP,iQ)*E(lT)
     & + 1.5D0*G1(iP,lT)*G1(lT,jQ)*G1(jP,iQ)*G1(kP,kQ)*E(lT)
     & + 1.5D0*G1(iP,lT)*G1(lT,kQ)*G1(kP,iQ)*G1(jP,jQ)*E(lT)
     & + 1.5D0*G1(jP,lT)*G1(lT,kQ)*G1(kP,jQ)*G1(iP,iQ)*E(lT)
c     // G(I,k)G(J,l)G(K,i)G(L,j)[3]
      CU4F3H=CU4F3H
     & + 1.5D0*G1(iP,kQ)*G1(jP,lT)*G1(kP,iQ)*G1(lT,jQ)*E(lT)
     & + 1.5D0*G1(iP,jQ)*G1(kP,lT)*G1(jP,iQ)*G1(lT,kQ)*E(lT)
     & + 1.5D0*G1(iP,lT)*G1(jP,kQ)*G1(lT,iQ)*G1(kP,jQ)*E(lT)
c     // G(I,j)G(J,k)G(K,l)G(L,i)[6]
      CU4F3H=CU4F3H
     & - 0.75D0*G1(iP,jQ)*G1(jP,kQ)*G1(kP,lT)*G1(lT,iQ)*E(lT)
     & - 0.75D0*G1(iP,kQ)*G1(kP,jQ)*G1(jP,lT)*G1(lT,iQ)*E(lT)
     & - 0.75D0*G1(iP,jQ)*G1(jP,lT)*G1(lT,kQ)*G1(kP,iQ)*E(lT)
     & - 0.75D0*G1(iP,kQ)*G1(kP,lT)*G1(lT,jQ)*G1(jP,iQ)*E(lT)
     & - 0.75D0*G1(iP,lT)*G1(lT,jQ)*G1(jP,kQ)*G1(kP,iQ)*E(lT)
     & - 0.75D0*G1(iP,lT)*G1(lT,kQ)*G1(kP,jQ)*G1(jP,iQ)*E(lT)
  900 CONTINUE
      Return
      End Function
