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
      Subroutine expandbas(Bas1,Nbas1,Bas2,Nbas2,Orb1,Orb2,occ1,eorb1,
     &                  indt1,occ2,eorb2,indt2)
C
C     This subroutine expands the MOs for a given symmetry
C     Orb1 are the input orbitals of dimension Nbas1*bas1 (file INPORB)
C     Orb2 are the output orbitals of dimension Nbas2*bas2 (file EXPORB)
C     Bas1 and Bas2 are the basis set specifications for the old and
C     new basis, respectively. They have dimensions Nbas1 and Nbas2.
C
C
      Implicit real*8 (a-h,o-z)
#include "Molcas.fh"
      Character*(LENIN4) Bas1(*),Bas2(*)
      Integer indt1(*),indt2(*)
      Dimension Orb1(*),Orb2(*),Izero(Nbas2),occ1(*),eorb1(*),
     &          occ2(*),eorb2(*)
C
C     Loop through the new basis labels and compare with the old.
C     If they are equal copy orbital coefficients
C     If not, add zeros until they fit again
C
      Nzero=0
      Ibas2=1
      If(Nbas1.eq.0) go to 200
      Ibas1=1
  100 Continue
      If(Bas2(Ibas2).eq.Bas1(Ibas1)) then
       lmo1=0
       lmo2=0
       Do imo=1,Nbas1
        Orb2(lmo2+Ibas2)=Orb1(lmo1+Ibas1)
        lmo1=lmo1+nBas1
        lmo2=lmo2+nBas2
       Enddo
       Ibas1=Ibas1+1
       Ibas2=Ibas2+1
      Else If(Bas2(Ibas2).ne.Bas1(Ibas1)) then
       Nzero=Nzero+1
       Izero(Nzero)=Ibas2
       lmo2=0
       Do imo=1,Nbas1
        Orb2(lmo2+Ibas2)=0.D0
        lmo2=lmo2+nBas2
       Enddo
       Ibas2=Ibas2+1
      Endif
      If(Ibas1.le.Nbas1) go to 100
  200 Continue
C
C     Add zeros at the end of each basis function
C
      If(Ibas2.le.Nbas2) then
       Do i=Ibas2,Nbas2
        Nzero=Nzero+1
        Izero(Nzero)=i
        lmo2=0
        Do imo=1,Nbas1
         Orb2(lmo2+i)=0.D0
         lmo2=lmo2+nBas2
        Enddo
       Enddo
      Endif
C
C     Add new basis functions at the end and expand occ, eorb, and indt
C
      If(nBas1.ne.0) then
       Do imo=1,nBas1
        occ2(imo)=occ1(imo)
        eorb2(imo)=eorb1(imo)
        indt2(imo)=indt1(imo)
       Enddo
      Endif

      If(Nbas1.lt.Nbas2) then
       Nzero=0
       Do imo=Nbas1+1,Nbas2
        Nzero=Nzero+1
        Do Ibas2=1,Nbas2
         Orb2(Nbas2*(imo-1)+Ibas2)=0.D0
        Enddo
        Orb2(Nbas2*(imo-1)+Izero(Nzero))=1.D0
        occ2(imo)=0.d0
        eorb2(imo)=0.d0
        indt2(imo)=6
       Enddo
      Endif

      Return
      End
