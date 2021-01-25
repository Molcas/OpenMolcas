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
      Subroutine FndBnd(IOut,IPrint,AlBond,ToAng,MxBond,
     +                  NAtoms,IAn,C,Nbond,IBond,
     +                  IBType,PBO,Re)
      Implicit Real*8(A-H,O-Z)
C
C     Generate connectivity based on bond distances alone.  The criteria
C     are contained in routine IPBO.
C     Cartesian coords. are in Angstroms
C
      Logical AlBond
C     Character AtSymb*2,AppNum*3
      Dimension IAn(NAtoms),C(3,NAtoms), NBond(NAtoms), Re(*),
     $  IBond(MxBond,NAtoms), IBType(MxBond,NAtoms), PBo(MxBond,NAtoms)
C     Data AtSymb, AppNum /'XX','XXX'/
 1000 Format(' Maximum number of bonds=',I3,' exceeded for atom',I4,'.')
C
      Do 10 I = 1, 12
      Do 11 J = 1, NAtoms
        IBond(I,J) = 0
        IBType(I,J) = 0
  11  Continue
  10  Continue
      BondOr = DBLE(0)
      Do 100 I = 1, NAtoms
        NBond(I) = 0
        Do 90 J = 1, NAtoms
          If(J.eq.I) goto 90
          RIJ = Sqrt((C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2+
     $               (C(3,I)-C(3,J))**2)
          IBondO = IPBO(ToAng,IAn(I),IAn(J),RIJ,BondOr)
          If(IBondO.gt.0.or.AlBond) then
            NBond(I) = NBond(I) + 1
            If(NBond(I).gt.MxBond) then
              Write(IOut,1000) MxBond, I
              Call Abend()
              endIf
            IBond(NBond(I),I) = J
            IBType(NBond(I),I) = IBondO
            PBO(NBond(I),I)=BondOr
            endIf
   90     Continue
  100   Continue
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(IPrint)
        Call Unused_real_array(Re)
      End If
      End
