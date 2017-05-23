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
      subroutine comp2ind(W,IDM,no)
      implicit none
      real*8 W,DD
      integer IDM,no, IJ,I,IJO,II,K,L,KL,LK, NF, J, JIO
      dimension W(IDM,*)
cmp!      include 'pvtrace_inc.f'
cmp!      include 'pvtrace_inc'
cmp!      if (LLtrace) then
cmp!         write(6,*) 'Entering COMP2IND'
cmp!         call xflush(6)
cmp!      endif
C
C   comprises the back square to triangle  (r,s,p,>t)>>> (r,s,pt)
C   remains in the same array!
C
C   Fixed DCOPY for no=2 to avoid overlap in input and output fields
C   with impredictable behavior on pwr4/ESSL.       PV, 12 aug 2004.
C
c      write(6,*)'comp2ind:',IDM,no
c   check
c      do I=1,no
c      do J=1,I
c      IJ=(I-1)*no+J
c      JI=(J-1)*no+I
c      DD1=DDOT_(IDM,W(1,IJ),1,W(1,IJ),1)
c      DD2=DDOT_(IDM,W(1,JI),1,W(1,JI),1)
c      if(abs(DD1-DD2).GT.1D-10)write(6,*)'warning'
c      write(6,*)dd1,dd2
c      enddo
c      enddo
      IJ=2
      do I=2,no
         IJO=(I-1)*no+1
C        CALL DCOPY_(IDM*I,W(1,IJO),1,W(1,IJ),1)
         if (no.eq.2) then
            CALL DCOPY_(IDM,W(1,3),1,W(1,2),1)
            CALL DCOPY_(IDM,W(1,4),1,W(1,3),1)
         else
            CALL DCOPY_(IDM*I,W(1,IJO),1,W(1,IJ),1)
         endif
         IJ=IJ+I
      enddo
cmp!      if (LLtrace) then
cmp!         write(6,*) 'Leaving COMP2IND'
cmp!         call xflush(6)
cmp!      endif
      return

      entry decomp2ind(W,IDM,no,NF)
cmp!      if (LLtrace) then
cmp!         write(6,*) 'Entering DECOMP2IND'
cmp!         call xflush(6)
cmp!      endif
C
C symmetrizes the upper index
C
C      write(6,*)'decomp2ind:',IDM,no,NF
      DO I=1,no
         II=I*(I+1)/2
         DO K=2,NF
            DO L=1,K-1
               KL=(K-1)*NF+L
               LK=(L-1)*NF+K
               DD=0.5d0*(W(KL,II)+W(LK,II))
               W(KL,II)=DD
               W(LK,II)=DD
            enddo
         enddo
      enddo
      IF(NO.GT.2)THEN
         DO I=NO,2,-1
            IJ=I*(I-1)/2+1
            IJO=(I-1)*no+1
            CALL DCOPY_(IDM*I,W(1,IJ),1,W(1,IJO),1)
         enddo
      ELSEIF(NO.EQ.2)THEN
         CALL DCOPY_(IDM,W(1,3),1,W(1,4),1)
         CALL DCOPY_(IDM,W(1,2),1,W(1,3),1)
      endif
      DO I=2,no
         DO J=1,I-1
            IJO=(I-1)*no+J
            JIO=(J-1)*no+I
            call transm(W(1,IJO),W(1,JIO),NF,NF)
         enddo
      enddo
c     stop
cmp!      if (LLtrace) then
cmp!         write(6,*) 'Leaving DECOMP2IND'
cmp!         call xflush(6)
cmp!      endif
      end
