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
C
C----------------------------------------------------------------------|
C
      subroutine dkh_wgene(n,ord,nst,ndk,ifodd,cdk,wr,rw,t1,t2,
     &                     e,rer,or,ro,info,s1,s2,t3,t4)
C
C Calculate U(Word)O/E_{nst}U^{\dag}(Word)
C
      implicit none
      integer n,ord,ndk,nst,info,m,i,j,k,L1,L2
      logical ifodd
      Real*8 cdk(ndk),t1(n,n),t2(n,n),e(n,n,ndk),rer(n,n,ndk),
     & or(n,n,ndk),ro(n,n,ndk),wr(n,n),rw(n,n),c
      Real*8 s1(n,n,*),s2(n,n,*),t3(n,n),t4(n,n)
C
      m=(ndk-nst)/ord+1
      if(m.le.1) return
      do L1=1,n
        do L2=1,n
          s1(L2,L1,1) = t1(L2,L1)
          s2(L2,L1,1) = t2(L2,L1)
        end do
      end do
      do i=1,m-1
        do L1=1,n
          do L2=1,n
            t1(L2,L1) = 0.d0
            t2(L2,L1) = 0.d0
          end do
        end do
        k=nst+ord*i
        call dkh_woprig(n,ifodd,ord,k-ord,wr,rw,s1(1,1,i),s2(1,1,i),
     &                    s1(1,1,i+1),s2(1,1,i+1),t3,t4 )
        info=info+2
        c = cdk(i) * (-1)**(i)
        do L1=1,n
          do L2=1,n
            t1(L2,L1) = t1(L2,L1) + s1(L2,L1,i+1) * c
            t2(L2,L1) = t2(L2,L1) + s2(L2,L1,i+1) * c
          end do
        end do
        do j=1,i
          call dkh_woplft(n,ifodd,ord,k-ord,wr,rw,s1(1,1,j),s2(1,1,j),
     &                      s1(1,1,j),s2(1,1,j),t3,t4 )
          info=info+2
          if(j.eq.1)then
            c = cdk(i)
          else
            c = cdk(i-j+1) * cdk(j-1) * (-1)**(j-1)
          end if
          do L1=1,n
            do L2=1,n
              t1(L2,L1) = t1(L2,L1) + s1(L2,L1,j) * c
              t2(L2,L1) = t2(L2,L1) + s2(L2,L1,j) * c
            end do
          end do
        end do
        ifodd=.not.ifodd
        if(ifodd)then
          do L1=1,n
            do L2=1,n
              or(L2,L1,k) = or(L2,L1,k) + t1(L2,L1)
              ro(L2,L1,k) = ro(L2,L1,k) + t2(L2,L1)
            end do
          end do
        else
          do L1=1,n
            do L2=1,n
              e(L2,L1,k)   = e(L2,L1,k)   + t1(L2,L1)
              rer(L2,L1,k) = rer(L2,L1,k) + t2(L2,L1)
            end do
          end do
        end if
      enddo
      end
