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
      subroutine mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *angintSO,angintOO,
     *Lfirst,Llast,Lblocks,
     *ncont1,ncont2,ncont3,
     *ncont4,
     *caseaSO,caseb1SO,caseb2SO,casecSO,
     *caseaOO,caseb1OO,caseb2OO,casecOO,
     *preroots,clebsch,dummy,bonn,breit,
     *sameorb)
      implicit real*8 (a-h,o-z)
cbs   subroutine for combining radial intgrls with angular
cbs   factors for the block with l1,l2,l3,l4,m1,m2,m3m,m4
cbs   this routine mkangLmin = make angular factors for the L- -part
cbs   includes both, spin-same and spin-other-orbit parts.
      real*8 LMdepang
      dimension
     *angintSO(ncont1,ncont2,ncont3,ncont4),
     *angintOO(ncont1,ncont2,ncont3,ncont4),
     *Lfirst(*),Llast(*),Lblocks(*),
cbs   all the arrays with the radial intgrls for
cbs   this combination of l-values
     *caseaSO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (2,0)   intgrls with alpha1*alpha3
     *caseb1SO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (0,0)   intgrls with alpha1
     *caseb2SO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (0,0)   intgrls with alpha3
     *casecSO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (-2,0)  intgrls with factor 1
     *caseaOO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (2,0)   intgrls with alpha1*alpha3
     *caseb1OO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (0,0)   intgrls with alpha1
     *caseb2OO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (0,0)   intgrls with alpha3
     *casecOO(ncont1*ncont2*ncont3*ncont4,*),
c    ! (-2,0)  intgrls with factor 1
     *preroots(2,0:Lmax),
c    ! some prefactors: sqrt( (l(+1))/(2l+1))
     *clebsch(3,2,-Lmax:Lmax,0:Lmax)
c    ! some clebsch gordans, that appear regulary
      dimension dummy(0:*)
      logical bonn,breiT,sameorb
      root2=sqrt(2.0d0)
      root2inv=1d0/root2
c     write(6,*) 'begin mkangL- ',
c    *l1,l2,l3,l4,m1,m2,m3,m4
cbs
      ncontall=ncont1*ncont2*ncont3*ncont4
cbs   cheater introduced to correct signs, because they were different from HERMIT
      if (mod(l1+l2+l3+l4,4).eq.2) then
      cheater=1d0
      else
      cheater=-1d0
      endiF
cbs   cleaning up
      if (bonn.or.breit.or.sameorb) then
      call dzero(angintSO,ncontall)
      else
      call dzero(angintSO,ncontall)
      call dzero(angintOO,ncontall)
      endif
cbs  starting with the same-orbit-contributions
cbs  first term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(2,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(2,2,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   4.0D0*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      Endif
      Lrun=Lrun+1
      enddo
      endif
      endif
cbs   second term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(2,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(2,2,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   4.0D0*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(3).gt.0) then
      M=m2-m4
              if (Lfirst(3).lt.Kfirst) then
              do L=Lfirst(3),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Kfirst=Lfirst(3)
              endif
              if (Llast(3).gt.Klast) then
              do L=Klast,Llast(3),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Klast=Llast(3)
              endif
      Lrun=1
      do L=Lfirst(3),Llast(3),2
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,-DBLE(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      endif
cbs   third term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(1,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(2,1,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   4.0D0*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(2).gt.0) then
      M=m2-m4
              if (Lfirst(2).lt.Kfirst) then
              do L=Lfirst(2),Kfirst,2
              dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,
     *                 m3,m4,Cheater)
              enddo
              Kfirst=Lfirst(2)
              endif
              if (Llast(2).gt.Klast) then
              do L=Klast,Llast(2),2
        dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Klast=Llast(2)
              endif
      Lrun=1
      do L=Lfirst(2),Llast(2),2
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,-DBLE(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      endif
cbs   fourth term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(1,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(2,1,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   4.0D0*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(2).gt.0) then
      M=m2-m4
              if (Lfirst(2).lt.Kfirst) then
              do L=Lfirst(2),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Kfirst=Lfirst(2)
              endif
              if (Llast(2).gt.Klast) then
              do L=Klast,Llast(2),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Klast=Llast(2)
              endif
      Lrun=1
      do L=Lfirst(2),Llast(2),2
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,-DBLE(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(3).gt.0) then
      M=m2-m4
              if (Lfirst(3).lt.Kfirst) then
              do L=Lfirst(3),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Kfirst=Lfirst(3)
              endif
              if (Llast(3).gt.Klast) then
              do L=Klast,Llast(3),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Klast=Llast(3)
              endif
      Lrun=1
      do L=Lfirst(3),Llast(3),2
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,-DBLE(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(4).gt.0) then
      M=m2-m4
              if (Lfirst(4).lt.Kfirst) then
              do L=Lfirst(4),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Kfirst=Lfirst(4)
              endif
              if (Llast(4).gt.Klast) then
              do L=Klast,Llast(4),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
              enddo
              Klast=Llast(4)
              endif
      Lrun=1
      do L=Lfirst(4),Llast(4),2
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,DBLE(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,DBLE(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   DBLE(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      endif
cbs  fifth term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(2,l3)*
     *clebsch(2,2,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   4.0D0*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      endif
cbs   sixth  term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(2,l3)*
     *clebsch(2,1,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,4.0D0*
     *   factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(3).gt.0) then
      M=m2-m4
              if (Lfirst(3).lt.Kfirst) then
              do L=Lfirst(3),Kfirst,2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Kfirst=Lfirst(3)
              endif
              if (Llast(3).gt.Klast) then
              do L=Klast,Llast(3),2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Klast=Llast(3)
              endif
      Lrun=1
      do L=Lfirst(3),Llast(3),2
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,-DBLE(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      endif
cbs   seventh term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(1,l3)*
     *clebsch(2,2,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   4.0D0*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(2).gt.0) then
      M=m2-m4
              if (Lfirst(2).lt.Kfirst) then
              do L=Lfirst(2),Kfirst,2
         dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Kfirst=Lfirst(2)
              endif
              if (Llast(2).gt.Klast) then
              do L=Klast,Llast(2),2
         dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Klast=Llast(2)
              endif
      Lrun=1
      do L=Lfirst(2),Llast(2),2
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,-DBLE(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      endif
cbs   eigth term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(1,l3)*
     *clebsch(2,1,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients
      if (Lblocks(1).gt.0) then
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,4.0D0*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,4.0D0*
     *   factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(2).gt.0) then
      M=m2-m4
              if (Lfirst(2).lt.Kfirst) then
              do L=Lfirst(2),Kfirst,2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Kfirst=Lfirst(2)
              endif
              if (Llast(2).gt.Klast) then
              do L=Klast,Llast(2),2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Klast=Llast(2)
              endif
      Lrun=1
      do L=Lfirst(2),Llast(2),2
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
        call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *  caseb1SO(1,Lrun),1,angintSO,1)
      else
        call daxpy_(ncontall,-DBLE(2+4*l3)*factor*dummy(L),
     *  caseb1SO(1,Lrun),1,angintSO,1)
        call daxpy_(ncontall,-DBLE(2+4*l3)*
     *factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(3).gt.0) then
      M=m2-m4
              if (Lfirst(3).lt.Kfirst) then
              do L=Lfirst(3),Kfirst,2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Kfirst=Lfirst(3)
              endif
              if (Llast(3).gt.Klast) then
              do L=Klast,Llast(3),2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Klast=Llast(3)
              endif
      Lrun=1
      do L=Lfirst(3),Llast(3),2
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,-DBLE(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,-DBLE(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      if (Lblocks(4).gt.0) then
      M=m2-m4
              if (Lfirst(4).lt.Kfirst) then
              do L=Lfirst(4),Kfirst,2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Kfirst=Lfirst(4)
              endif
              if (Llast(4).gt.Klast) then
              do L=Klast,Llast(4),2
         dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
              enddo
              Klast=Llast(4)
              endif
      Lrun=1
      do L=Lfirst(4),Llast(4),2
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then
         call daxpy_(ncontall,DBLE(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else
         call daxpy_(ncontall,DBLE(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy_(ncontall,
     *   DBLE(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif
      endif
      Lrun=Lrun+1
      enddo
      endif
      endif
      return
      end
