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
      Subroutine priwel(k,alfa,beta,r0,a,gri,nz,isum,grin)
      Implicit Real*8(A-H,O-Z)
#include "real.fh"
#include "welcom.fh"
      Real*8 gri(nz,isum), grin(nz,0:k,k/2+1,k/4+1), alfa(nz), a(nz)
      Integer iv(kmax)
*
*     iQ = 1
      Call qEnter('PriWel')
      Call binte(k,alfa,beta,r0,a,grin,nz)
*     Call RecPrt(' In PriWel: Grin',' ',Grin,nz,(k+1)*(k/2+1)*(k/4+1))
*
c.....distribute the integrals into gri
*
      indst=1
      isum=ipot3(k+1)
      Do 115 iz=1,nz
         gri(iz,1)=grin(iz,0,1,1)
 115  Continue
      if (k.eq.0) Go To 99
      Do 10 i=1,k
         ipot3i=ipot3(i)
         call dcopy_(iPot3i*nz,Zero,0,gri(1,indst+1),1)
         Do 11 j=1,ipot3i
            jj=j
            Do 12 l=i,1,-1
               idiv=ipot3(l-1)
               ixyz=(jj-1)/idiv+1
               iv(l)=ixyz
               jj=jj-(ixyz-1)*idiv
12          Continue
*
c.....the potency vector for this integral is now ready
c .....now analyze it
*
            ix=0
            iy=0
            iz=0
            Do 13 l=1,i
               if(iv(l).eq.1)ix=ix+1
               if(iv(l).eq.2)iy=iy+1
               if(iv(l).eq.3)iz=iz+1
 13         Continue
            ix2=(ix/2)*2
            if(ix2.ne.ix) go to 11
            iy2=(iy/2)*2
            if(iy2.ne.iy)go to 11
            ixs=max(ix,iy)
            iys=min(ix,iy)
            ixys=(ixs+iys)/2+1
            iys=iys/2+1
            Do 14 mz=1,nz
               gri(mz,j+indst)=grin(mz,i,ixys,iys)
 14         Continue
 11      Continue
         indst=indst+ipot3i
 10   Continue
 99   Continue
*     Call RecPrt(' In PriWel:gri',' ',gri,nz,isum)
*
      Call qExit('PriWel')
      Return
      End
