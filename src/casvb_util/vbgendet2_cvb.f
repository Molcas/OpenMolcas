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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine vbgendet2_cvb(
     >  iapr,ixapr,ibpr,ixbpr,
     >  iconfs,idetvb,
     >  nconf,nconfion,
     >  nda,ndb,ndetvb,nel,
     >  noe,nalf,nbet,norb,
     >  idetavb,idetbvb,iwrk1,iwrk2)
      implicit real*8 (a-h,o-w,y-z),integer(x)
      dimension iapr(ndetvb),ixapr(nda+1),ibpr(ndetvb),ixbpr(ndb+1)
      dimension iconfs(noe,nconf),idetvb(ndetvb)
      dimension nconfion(0:nel)
      dimension idetavb(ndetvb),idetbvb(ndetvb)
      dimension iwrk1(ndetvb),iwrk2(ndetvb)
      logical debug
      data debug/.false./

      if(debug)then
        write(6,*)' Generate determinant information :'
        write(6,*)' ----------------------------------'
      endif

c  vbgenabdet gives all VB alpha and beta strings, to use in CASSCF
c  space we sort in A/B strings to get increasing order :
      call vbgenabdet_cvb(idetavb,idetbvb,
     >  iconfs,nconf,nconfion,
     >  ndetvb,nel,noe,
     >  nalf,nbet,norb)

      call sortindxi_cvb(ndetvb,idetbvb,idetvb)
      do 100 i=1,ndetvb
      iwrk1(i)=idetbvb(idetvb(i))
      ibpr(i)=idetavb(idetvb(i))
100   continue
      ixbpr(1)=1
      do 200 i=1,ndb
      do 300 j=ixbpr(i),ndetvb
      if(iwrk1(j).ne.i)goto 400
300   continue
      j=ndetvb+1
400   ixbpr(i+1)=j
200   continue
      do 500 i=1,ndb
      call sortindxi_cvb(ixbpr(i+1)-ixbpr(i),ibpr(ixbpr(i)),iwrk2)
      do 600 j=1,ixbpr(i+1)-ixbpr(i)
      iwrk1(j)=ibpr(iwrk2(j)+ixbpr(i)-1)
600   continue
      call imove_cvb(iwrk1,ibpr(ixbpr(i)),ixbpr(i+1)-ixbpr(i))
500   continue
      if(debug)then
        write(6,*)' ixbpr='
        write(6,'(10i6)')ixbpr
        write(6,*)' ibpr='
        write(6,'(10i6)')ibpr
      endif

      call sortindxi_cvb(ndetvb,idetavb,idetvb)
      do 700 i=1,ndetvb
      iwrk1(i)=idetavb(idetvb(i))
      iapr(i)=idetbvb(idetvb(i))
700   continue
      ixapr(1)=1
      do 800 i=1,nda
      do 900 j=ixapr(i),ndetvb
      if(iwrk1(j).ne.i)goto 1000
900   continue
      j=ndetvb+1
1000  ixapr(i+1)=j
800   continue
      do 1100 i=1,nda
      call sortindxi_cvb(ixapr(i+1)-ixapr(i),iapr(ixapr(i)),iwrk2)
      do 1200 j=1,ixapr(i+1)-ixapr(i)
      iwrk1(j)=iapr(iwrk2(j)+ixapr(i)-1)
1200  continue
      call imove_cvb(iwrk1,iapr(ixapr(i)),ixapr(i+1)-ixapr(i))
      do 1300 j=1,ixapr(i+1)-ixapr(i)
      iwrk1(j)=idetvb(iwrk2(j)+ixapr(i)-1)
1300  continue
      call imove_cvb(iwrk1,idetvb(ixapr(i)),ixapr(i+1)-ixapr(i))
1100  continue
      if(debug)then
        write(6,*)' ixapr='
        write(6,'(10i6)')ixapr
        write(6,*)' iapr='
        write(6,'(10i6)')iapr
      endif
      do 1400 i=1,ndetvb
      idetavb(idetvb(i))=i
1400  continue
      call imove_cvb(idetavb,idetvb,ndetvb)
      if(debug)then
        write(6,*)' idetvb='
        write(6,'(10i6)')idetvb
      endif
      return
      end
