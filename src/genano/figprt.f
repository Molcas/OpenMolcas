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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
* This routine prints the body part to the postscript figure of        *
* occupation numbers.                                                  *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine FigPrt(lu,Title,nLqn,nOrb,Occ)
      Implicit Real*8 (a-h,o-z)
      Dimension nOrb(0:nLqn),Occ(*)
      Dimension iOff(0:25)
      Character*(*) Title
      Character*1 type(0:5)
      Data type/'s','p','d','f','g','h'/
      Write(lu,'(a)') 'HF setfont'
      Write(lu,'(a)') '2.25 X -7.0 Y moveto'
      k0 =0
      k1 =0
      Do 10 i=len(Title),1,-1
         If(Title(i:i).ne.' ') k0=i
10    Continue
      Do 20 i=1,len(Title)
         If(Title(i:i).ne.' ') k1=i
20    Continue
      Write(lu,'(3a)') '(',Title(k0:k1),') CenterLine'
      nOff=0
      mxNqn=0
      Do 100 i=0,nLqn
         if(nOrb(i).gt.0) mxNqn=max(mxNqn,nOrb(i)+i)
         iOff(i)=nOff
         nOff=nOff+nOrb(i)
100   Continue
      Do 200 n=1,mxNqn
         Write(lu,'(a,i2)') '%--- shell n=',n
         lqn0=nLqn
         lqn1=0
         Do 210 i=0,min(n-1,nLqn)
            If(n-i.le.nOrb(i)) Then
               If(i.gt.lqn1) lqn1=i
               If(i.lt.lqn0) lqn0=i
               O=log10(Occ(iOff(i)+n-i))
               Write(lu,'(i2,a,f7.4,a)') i,' X ',O,' Y Draw'
            End If
210      Continue
         if(n.lt.10) Then
            Write(lu,'(a,i1,a,a,i1,a,a)')
     &         ' (',n,type(lqn0),'\261',n,type(lqn1),') Label'
         Else
            Write(lu,'(a,i2,a,a,i2,a,a)')
     &         ' (',n,type(lqn0),'\261',n,type(lqn1),') Label'
         End If
         Do 220 i=1,min(n-1,nLqn)
            If( (n-i.le.nOrb(i)).and.(n-i+1.le.nOrb(i-1)) ) Then
               Oa=log10(Occ(iOff(i-1)+n-i+1))
               Ob=log10(Occ(iOff(i)+n-i))
               Write(lu,'(i2,a,f7.4,a,i2,a,f7.4,a)')
     &            i-1,' X ',Oa,' Y ',i,' X ',Ob,' Y Connect'
            End If
220      Continue
200   Continue
      Return
      End
