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
* Copyright (C) Ben Swerts                                             *
*               2016, Liviu Ungur                                      *
************************************************************************
      Subroutine FragpCont(F1,mi,mK,ma,mC,F2,mL,mj,mD,mb,W,Final,Factor)
************************************************************************
*                                                                      *
* Object: Specialized contraction of 3 4D matrices.                    *
*         It is very inefficient and should be rewritten using BLAS    *
*                                                                      *
* Called from: FragPInt and FragPGrd                                   *
*                                                                      *
*     Author: Ben Swerts                                               *
*   Modified: Liviu Ungur                                              *
************************************************************************
      Implicit None
      Integer mi,mK,ma,mC
      Integer mL,mj,mD,mb
      Real*8     F1( mi, mK, ma, mC )
      Real*8     F2( mL, mj, mD, mb )
      Real*8      W( mK, mC, mL, mD )
      Real*8  Final( mi, mj, ma, mb )
      Real*8  Factor
c local variables
      Integer ib,ia,ij,ii,iC,iK,iD,iL
      Real*8  xt
      Logical DBG
cc
      Integer j1,j2,j12
      Real*8 xt2, WW1(mK*mC*mL*mD),F12(mK*mC*mL*mD)
      Real*8 ddot_, dnrm2_
      External ddot_, dnrm2_

cc
      DBG=.false.
cc
      if(DBG) write(6,'(A)') 'Enter FragpCont'
      if(DBG) write(6,'(A,i5)') 'FragpCont:  mi=',mi
      if(DBG) write(6,'(A,i5)') 'FragpCont:  mK=',mK
      if(DBG) write(6,'(A,i5)') 'FragpCont:  ma=',ma
      if(DBG) write(6,'(A,i5)') 'FragpCont:  mC=',mC
      if(DBG) write(6,'(A,i5)') 'FragpCont:  mL=',mL
      if(DBG) write(6,'(A,i5)') 'FragpCont:  mj=',mj
      if(DBG) write(6,'(A,i5)') 'FragpCont:  mD=',mD
      if(DBG) write(6,'(A,i5)') 'FragpCont:  mb=',mb
      if(DBG) write(6,'(A,F20.14)') 'FragpCont:  Factor=',Factor
c
      if(DBG) write(6,'(A,i5)') 'FragpCont: F1(mi,mK,ma,mC):'
      if(DBG)  then
      do ii=1,mi
        do iK=1,mK
          do ia=1,ma
            do iC=1,mC
            write(6,'(4(a,i3),a,e22.14)')
     & 'FragpCont: F1(',ii,',',iK,',',ia,',',iC,')=',F1(ii,iK,ia,iC)
            enddo
          enddo
        enddo
      enddo
c
      write(6,'(A,i5)') 'FragpCont: F2(mL,mj,mD,mb):'
      do iL=1,mL
        do ij=1,mj
          do iD=1,mD
            do ib=1,mb
            write(6,'(4(a,i3),a,e22.14)')
     & 'FragpCont: F2(',iL,',',ij,',',iD,',',ib,')=',F2(iL,ij,iD,ib)
            enddo
          enddo
        enddo
      enddo
c
      write(6,'(A,i5)') 'FragpCont: W(mK,mC,mL,mD):'
      do iK=1,mK
        do iC=1,mC
          do iL=1,mL
            do iD=1,mD
            write(6,'(4(a,i3),a,e22.14)')
     & 'FragpCont:  W(',iK,',',iC,',',iL,',',iD,')=',W(iK,iC,iL,iD)
            enddo
          enddo
        enddo
      enddo
      endif !DBG

c      print *, 'norm of FINAL in fragpcont',
c     &          dnrm2_(mi*mj*ma*mb, Final(1:mi,1:mj,1:ma,1:mb), 1 )
      ! T = F1*W
c              Call DGEMM_('T','N',
c     &                    iBas*nElem(la)*nAlpha,iSize,nElem(iAng),
c     &                    1.0d0,Array(ipTmp),nElem(iAng),
c     &                    RSph(ipSph(iAng)),nElem(iAng),
c     &                    0.0d0,Array(ipF1),nAlpha*iBas*nElem(la))
      ! Final=T*F2
c              Call DGEMM_('T','N',
c     &                    iBas*nElem(la)*nAlpha,iSize,nElem(iAng),
c     &                    1.0d0,Array(ipTmp),nElem(iAng),
c     &                    RSph(ipSph(iAng)),nElem(iAng),
c     &                    0.0d0,Array(ipF1),nAlpha*iBas*nElem(la))
c      print *,'factor=', factor
      do ib = 1, mb
        do ia = 1, ma
          do ij = 1, mj
            do ii = 1, mi
            xt = 0.d0
            xt2= 0.d0
            WW1= 0.d0
            F12= 0.d0


              j1=0
              j12=0
              do iC = 1, mC
                do iK = 1, mK
                j1=j1+1

                  j2=0
                  do iD = 1, mD
                    do iL = 1, mL
                   j2=j2+1
                   j12=j12+1

                 WW1(j12)=W(iK,iC,iL,iD)
                 F12(j12)=F1(ii,iK,ia,iC)*F2(iL,ij,iD,ib)

            xt = xt + F1(ii,iK,ia,iC) * W(iK,iC,iL,iD) * F2(iL,ij,iD,ib)
                    enddo
                  enddo
                enddo
              enddo

              xt2=xt2+ddot_(mC*mK*mD*mL, F12, 1, WW1, 1)

c              print *, 'factor, xt, xt2, diff,:',factor,xt,xt2,xt-xt2

c              TMP
c            call dgemm_('N',  m, n,
c     &                 1.d0, F1(ii,:,ia,:), lda,
c     &                       W(:,:,iL,iD), incx,
c     &                  beta, TMP, incy )
c
c            call daxpy_(n, a, x, incx, y, incy)


            Final(ii,ij,ia,ib) = Final(ii,ij,ia,ib) + Factor * xt

            enddo !ii
          enddo !ij
        enddo !ia
      enddo !ib

      Return
      End
