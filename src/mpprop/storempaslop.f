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
      Subroutine StoreMpAsLop(nAtoms,ip_ANr,nB,ipT,ipTi
     &                       ,ipMP,lMax,ip_EC)
      Implicit Real*8 (a-h,o-z)

#include "MpParam.fh"
#include "Address.fh"
#include "WrkSpc.fh"
#include "MolProp.fh"

*
*-- Let's fix the ip_ANr.
*
      Call Allocate_iWork(ip_ANr,nAtoms)
      Call Get_iArray('LP_A',iWork(ip_ANr),nAtoms)

*
*-- Let's fix the uber-simple T and T(-1).
*
      Call GetMem('T','Allo','Real',ipT,nB**2)
      Call GetMem('Tinv','Allo','Real',ipTi,nB**2)
      kaunter=0
      Do i=1,nB
        Do j=1,nB
          Work(ipT+kaunter)=0.0d0
          Work(ipTi+kaunter)=0.0d0
          If(i.eq.j) Work(ipT+kaunter)=1.0d0
          If(i.eq.j) Work(ipTi+kaunter)=1.0d0
          kaunter=kaunter+1
        Enddo
      Enddo

*
*-- Let's fix the expansion centres. Cor resides in MolProp.fh.
*
      Call GetMem('ExpCent','Allo','Real',ip_EC,3*nAtoms*(nAtoms+1)/2)
      kaunter=0
      Do i=1,nAtoms
        Do j=1,i
          Work(ip_EC+kaunter*3+0)=Cor(1,i,j)
          Work(ip_EC+kaunter*3+1)=Cor(2,i,j)
          Work(ip_EC+kaunter*3+2)=Cor(3,i,j)
          kaunter=kaunter+1
        Enddo
      Enddo

*
*-- Let's fix the multipole moments. Unlike LoProp, MpProp has here
*   included the nuclei contribution, which we have to remove pronto
*   to be compatible.
*
      nSize1=nAtoms*(nAtoms+1)/2
      nSize2=(lMax*(lMax**2+6*lMax+11)+6)/6
      Call GetMem('MultMom','Allo','Real',ipMP,nSize1*nSize2)
      iMu=-1
      Do l=0,lMax
        kompost=0
        Do ix=l,0,-1
          Do iy=l-ix,0,-1
            kompost=kompost+1
            iMu=iMu+1
            iAtK=0
            Do iAt1=1,nAtoms
              Do iAt2=1,iAt1
                Work(ipMP+iAtK+nSize1*iMu)=
     &          Work(iAtBoMltPlAd(l)+nSize1*(kompost-1)+iAtK)
                iAtK=iAtK+1
              Enddo
              If(l.eq.0) then
                Work(ipMP+iAtK-1+nSize1*iMu)=
     &           Work(ipMP+iAtK-1+nSize1*iMu)-dble(iWork(ip_ANr+iAt1-1))
              Endif
            Enddo
          Enddo
        Enddo
      Enddo

      Return
      End
