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
      Subroutine HaveWeConv(iCNum,iCStart,iQ_Atoms,Indma,iDT,FFp,xyzMyI
     &,Egun,Energy,NVarv,JaNej,Haveri)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"

      Dimension iDT(3)
      Dimension FFp(npol*npart,3),xyzMyI(3)
      Logical JaNej,Haveri

*----------------------------------------------------------------------*
* With the new and the old induced dipoles, check if we have converged.*
* We also have energy check.                                           *
*----------------------------------------------------------------------*
      JaNej=.true.
      Haveri=.false.
      Diffab=0
      xyzMyi(1)=0
      xyzMyi(2)=0
      xyzMyi(3)=0
      Do 821, i=1+(nPol*iCnum),IndMa
        k=i-((i-1)/nPol)*nPol
        Do 822, l=1,3
          Dtil=FFp(i,l)*Pol(k)
          Diff=Abs(Work(iDT(l)+i-1)-Dtil)
          If(Diff.gt.Diffab) Diffab=Diff
          xyzMyi(l)=xyzMyi(l)+Dtil
          Work(iDT(l)+i-1)=Dtil  !This is the quantities that has
                !changed during the iteration and that through FFp
                !includes the effect of the polarization of the
                !QM-molecule. It enters the iteration above, unless
                !we have converged.
822     Continue
821   Continue
      Egtest=Egun-Energy
      Egun=Energy
      If(nVarv.ge.itMax) then  !itMax is from input or default.
        Write(6,*)
        Write(6,*)'  No convergence for the induced dipoles.'
        Write(6,*)'  Difference remaining after ',nVarv,' iterations: '
     &,Diffab
        Haveri=.true.
        iPrint=10
        Do 842, j=icstart,npart*ncent,ncent
          distmin=1000.0d0
          kmin=0
          imin=0
          Do 841, i=1,iq_atoms
            Do 843, k=0,ncent-1
              dist=sqrt((cordst(i,1)-cordst(j+k,1))**2
     &            +(cordst(i,2)-cordst(j+k,2))**2
     &            +(cordst(i,3)-cordst(j+k,3))**2)
              if(dist.lt.distmin) then
                distmin=dist
                imin=i
                kmin=k
              endif
843         Continue
841       Continue
         Write(6,*)'solv.',j,'iq_atom',imin,'center',kmin+1
     &            ,'dist',distmin
842     Continue
        Write(6,*)
        GoTo 9898
      Endif
      If(abs(egtest).gt.Enelim) JaNej=.false.
      If(Diffab.gt.PolLim) JaNej=.false.

9898  Continue

      Return
      End
