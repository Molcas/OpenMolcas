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
      Subroutine Local_XHole(ipXHole2,dMolExpec,nAtoms,nBas1,nTemp
     &                      ,iCenter,Ttot,Ttot_Inv,Coor
     &                      ,nij,EC,iANr,Bond_Threshold,iPrint
     &                      ,ipXHLoc2)
      Implicit Real*8 (a-h,o-z)

#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Dimension Coor(3,nAtoms), Sq_Temp(nTemp), Ttot_Inv(nTemp)
      Dimension Temp(nTemp), A(3), B(3), d2Loc(nij), EC(3,nij)
      Dimension Ttot(nTemp)
      Dimension iCenter(nBas1), iANr(nAtoms)
      Logical Found
      Real*8, Allocatable:: Dens(:)

*
*---- Binomal stuff
*
      Call Set_Binom()

*
*---- Transform the density matrix to the LoProp basis
*
      Call Qpg_dArray('D1ao',Found,nDenno)
      If (Found .and. nDenno/=0) Then
         Call mma_Allocate(Dens,nDenno,Label='Dens')
      Else
         Write (6,*) 'Local XHole: D1ao not found.'
         Call Abend()
      End If
      Call Get_D1ao(Dens,nDenno)
      Call DSq(Dens,Sq_Temp,1,nBas1,nBas1)
      Call mma_deallocate(Dens)
*
      Call DGEMM_('N','T',
     &            nBas1,nBas1,nBas1,
     &            1.0d0,Sq_Temp,nBas1,
     &            Ttot_Inv,nBas1,
     &            0.0d0,Temp,nBas1)
      Call DGEMM_('N','N',
     &            nBas1,nBas1,nBas1,
     &            1.0d0,Ttot_Inv,nBas1,
     &            Temp,nBas1,
     &            0.0d0,Sq_Temp,nBas1)

*
*---- Transform the exchange-hole stuff to LoProp basis
*
      Call Transmu(Work(ipXHole2),nBas1,Ttot,Temp)

*
*---- Now localize
*
      Do iAtom = 1, nAtoms
        call dcopy_(3,Coor(1,iAtom),1,A,1)
        Do jAtom = 1, iAtom
          call dcopy_(3,Coor(1,jAtom),1,B,1)
*
*-------- Sum up contibutions to the domain ij
*
          Acc=Zero
          iOffO=ipXHole2-1
          Do j = 1, nBas1
            Do i = 1, nBas1
              If ((iCenter(i).eq.iAtom .and.
     &             iCenter(j).eq.jAtom) .or.
     &            (iCenter(i).eq.jAtom .and.
     &             iCenter(j).eq.iAtom)) Then
                ij = (j-1)*nBas1+i
                Acc = Acc + Sq_Temp(ij)*Work(ij+iOffO)
              End If
            End Do
          End Do
          ij=iAtom*(iAtom-1)/2+jAtom
          d2Loc(ij)=Acc
        End Do   ! jAtom
      End Do   ! iAtom

*
* Distributes the contributions from the bonds that doesn't fulfill the requirement
* Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
* two atoms involved in the bond.
*
      Call Move_XHole(d2Loc,EC,nAtoms,nij,iANr,Bond_Threshold)
      call dcopy_(nij,d2Loc,1,Work(ipXHLoc2),1)

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(dMolExpec)
         Call Unused_integer(iPrint)
      End If
      End
