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
      Subroutine ReaInd(iGP,iDT,iDistIm,iCNum,indma,ncparm,Sum1,s90um)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "WrkSpc.fh"

      Dimension iGP(3),iDT(3)

      Sum1=0
*      irekn=0
*      xled=0
*      yled=0
*      zled=0
      Do 901, i=1+(nPol*iCnum),indma
        Do 902, j=1,3  !The energy of the induced dipoles (iDT) in
                       !the field from the real charges. Yes, this
                       !is the polarization energy in a system
                       !without permanent dipoles, see good old
                       !Bottcher, eq. (3.129). Observe also that
                       !we here include effects of the reaction
                       !field on the induced dipoles, see the
                       !polarization loop.
          Sum1=Sum1+Work(iGP(j)+i-1)*Work(iDT(j)+i-1)
902     Continue
* IF WE WISH TO MONITOR THE INDUCED DIPOLES, UNCOMMENT THIS, AND THE
* COMMENTED THING ABOVE.
*        irekn=irekn+1
*        xled=xled+Work(iDt(1)+i-1)
*        yled=yled+Work(iDt(2)+i-1)
*        zled=zled+Work(iDt(3)+i-1)
*        if(irekn.eq.3) then
*          irekn=0
*          TOT=sqrt(xled**2+yled**2+zled**2)
*          write(6,*)'HHH',TOT
*          xled=0
*          yled=0
*          zled=0
*        endif
901   Continue
      Sum1=Sum1*0.5
      S90um=0
      Do 911, i=iCnum+1,nPart !Energy of charge distribution
        Do 912, j=1,nPol    !in the reaction field to
          Q1=Qimp((i-1)*nPol+j) !the induced dipoles.
          D1x=Dim((i-1)*nPol+j,1) !Once more, see Bottcher eq.
          D1y=Dim((i-1)*nPol+j,2) !(4.69): we are computing the
          D1z=Dim((i-1)*nPol+j,3) !product between charges and the
          x=CordIm((i-1)*nCent+j,1) !potential connected with the
          y=CordIm((i-1)*nCent+j,2)  !reaction field.
          z=CordIm((i-1)*nCent+j,3)
          Inc=ncparm*nCent*(i-(iCnum+1))+(j-1)*ncparm
          Do 913, l=nCent-nCha+1,nCent
            Inc2=Inc+l
            Q2=Qsta(l-nCent+nCha)
            Do 914, k=iCnum+1,nPart
              X1=(X-Cordst(l+(k-1)*nCent,1))*D1x
              X1=(Y-Cordst(l+(k-1)*nCent,2))*D1y+X1
              X1=(Z-Cordst(l+(k-1)*nCent,3))*D1z+X1
         !Change sign on Q2 since we are in the backwards land, while
         !Q1 and X1 already are backward.
         S90um=S90um-(Q1+X1*Work(iDistIm-1+inc2+(k-(iCnum+1))*nCent)**2)
     &        *Q2*Work(iDistIm-1+inc2+(k-(iCnum+1))*nCent)
914         Continue
913       Continue
912     Continue
911   Continue

      Return
      End
