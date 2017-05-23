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
      Subroutine SolveA(AlfMat,AlfMatI,dLambda,dMullig,lMax
     &                 ,ARaw,BRaw,dA,iPrint,AboveMul,ddUpper,ddLower)
      Implicit Real*8 (a-h,o-z)

      Dimension AlfMat(4),AlfMatI(4),ARaw(2,2),BRaw(2),Beta(2),dA(2)
      Dimension dtA(2),dMullig((lMax*(lMax**2+6*lMax+11)+6)/6)

      Logical AboveMul(2)
      Logical Yeps(2)

*
*-- Which elements can be non-zero? Too low magnitude of multipoles
*   and they are screened.
*
      nDim=0
      Do iMull=1,2
        If(AboveMul(iMull)) then
          Yeps(iMull)=.true.
          nDim=nDim+1
          Beta(nDim)=BRaw(iMull)
        Else
          Yeps(iMull)=.false.
        Endif
      Enddo
      If(iPrint.ge.10) then
        Call RecPrt('Beta',' ',Beta,nDim,1)
      Endif

*
*-- Shuffle around, and create a matrix for exponents with non-zero
*   factors.
*
      kaunt=0
      Do i=1,2
        Do j=1,2
          If(Yeps(i).and.Yeps(j)) then
            kaunt=kaunt+1
            If(i.eq.j) then
              AlfMat(kaunt)=ARaw(max(i,j),min(i,j))*(1.0d0+dLambda)
            Else
              AlfMat(kaunt)=ARaw(max(i,j),min(i,j))
            Endif
          Endif
        Enddo
      Enddo

*
*-- Invert and solve.
*
      Call MInv(AlfMat,AlfMatI,Ising,Det,nDim)
      call dcopy_(nDim,0.0d0,0,dtA,1)
      Call dGeMV_('N',nDim,nDim,1.0d0,AlfMatI,nDim,Beta,1,0.0d0,dtA,1)

*
*-- Optional printing.
*
      If(iPrint.ge.10) then
        Call RecPrt('Alfa',' ',AlfMat,nDim,nDim)
        Call RecPrt('InverseA',' ',AlfMatI,nDim,nDim)
        Call RecPrt('deltatA',' ',dtA,nDim,1)
      Endif

*
*-- Damp large steps since such steps can take the optimization
*   too far away and put it in a region with very small derivatives,
*   and there things turns into baloney.
*
      If(dtA(1).lt.ddLower) dtA(1)=ddLower
      If(dtA(2).lt.ddLower) dtA(2)=ddLower
      If(dtA(1).gt.ddUpper) dtA(1)=ddUpper
      If(dtA(2).gt.ddUpper) dtA(2)=ddUpper

*
*-- Extend to full dimension.
*
      kaunt=0
      Do i=1,2
        If(Yeps(i)) then
          kaunt=kaunt+1
          dA(i)=dtA(kaunt)
        Else
          dA(i)=0.0d0
        Endif
      Enddo

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(dMullig)
      End
