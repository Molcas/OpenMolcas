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
      subroutine ass3b(D1,D2,D3,PAO,tmp1_,tmp2_,tmp3_,nt,nrys)
      implicit real*8 (a-h,o-z)
      dimension D1(nrys,nt),d2(nrys,nt),D3(nrys,nt),PAO(nt)
c
      If (nRys.eq.1) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + PAO(iT) * D1(1,iT)
        tmp2 = tmp2 + PAO(iT) * D2(1,iT)
        tmp3 = tmp3 + PAO(iT) * D3(1,iT)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return

      Else If (nRys.eq.2) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
c
      Else If (nRys.eq.3) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)
     >              +  D1(3,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)
     >              +  D2(3,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)
     >              +  D3(3,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
c
      Else If (nRys.eq.4) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)
     >              +  D1(3,iT)
     >              +  D1(4,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)
     >              +  D2(3,iT)
     >              +  D2(4,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)
     >              +  D3(3,iT)
     >              +  D3(4,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return

      Else If (nRys.eq.5) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)
     >              +  D1(3,iT)
     >              +  D1(4,iT)
     >              +  D1(5,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)
     >              +  D2(3,iT)
     >              +  D2(4,iT)
     >              +  D2(5,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)
     >              +  D3(3,iT)
     >              +  D3(4,iT)
     >              +  D3(5,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
c
      Else If (nRys.eq.6) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)
     >              +  D1(3,iT)
     >              +  D1(4,iT)
     >              +  D1(5,iT)
     >              +  D1(6,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)
     >              +  D2(3,iT)
     >              +  D2(4,iT)
     >              +  D2(5,iT)
     >              +  D2(6,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)
     >              +  D3(3,iT)
     >              +  D3(4,iT)
     >              +  D3(5,iT)
     >              +  D3(6,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
c
      Else If (nRys.eq.7) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)
     >              +  D1(3,iT)
     >              +  D1(4,iT)
     >              +  D1(5,iT)
     >              +  D1(6,iT)
     >              +  D1(7,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)
     >              +  D2(3,iT)
     >              +  D2(4,iT)
     >              +  D2(5,iT)
     >              +  D2(6,iT)
     >              +  D2(7,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)
     >              +  D3(3,iT)
     >              +  D3(4,iT)
     >              +  D3(5,iT)
     >              +  D3(6,iT)
     >              +  D3(7,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
c
      Else If (nRys.eq.8) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)
     >              +  D1(3,iT)
     >              +  D1(4,iT)
     >              +  D1(5,iT)
     >              +  D1(6,iT)
     >              +  D1(7,iT)
     >              +  D1(8,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)
     >              +  D2(3,iT)
     >              +  D2(4,iT)
     >              +  D2(5,iT)
     >              +  D2(6,iT)
     >              +  D2(7,iT)
     >              +  D2(8,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)
     >              +  D3(3,iT)
     >              +  D3(4,iT)
     >              +  D3(5,iT)
     >              +  D3(6,iT)
     >              +  D3(7,iT)
     >              +  D3(8,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
c
      Else If (nRys.eq.9) Then
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      do it=1,nt
        tmp1 = tmp1 + (D1(1,iT)
     >              +  D1(2,iT)
     >              +  D1(3,iT)
     >              +  D1(4,iT)
     >              +  D1(5,iT)
     >              +  D1(6,iT)
     >              +  D1(7,iT)
     >              +  D1(8,iT)
     >              +  D1(9,iT)) * PAO(It)
        tmp2 = tmp2 + (D2(1,iT)
     >              +  D2(2,iT)
     >              +  D2(3,iT)
     >              +  D2(4,iT)
     >              +  D2(5,iT)
     >              +  D2(6,iT)
     >              +  D2(7,iT)
     >              +  D2(8,iT)
     >              +  D2(9,iT)) * PAO(It)
        tmp3 = tmp3 + (D3(1,iT)
     >              +  D3(2,iT)
     >              +  D3(3,iT)
     >              +  D3(4,iT)
     >              +  D3(5,iT)
     >              +  D3(6,iT)
     >              +  D3(7,iT)
     >              +  D3(8,iT)
     >              +  D3(9,iT)) * PAO(It)
      enddo
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
c
      Else
      tmp1=0.0D0
      tmp2=0.0D0
      tmp3=0.0D0
      Do iT = 1, nT
        Do iRys = 1, nRys
          tmp1 = tmp1 + PAO(It)*D1(iRys,iT)
          tmp2 = tmp2 + PAO(It)*D2(iRys,iT)
          tmp3 = tmp3 + PAO(It)*D3(iRys,iT)
        End Do
      End Do
      tmp1_ = tmp1_ + tmp1
      tmp2_ = tmp2_ + tmp2
      tmp3_ = tmp3_ + tmp3
      return
      End If
      end
