!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!----------------------------------------------------------------------*
! A subroutine that computes those darn f-factors. They are definied   *
! in equation (2.4) in the article cited above. As can be seen from    *
! that equation, the computation of the f-factors is actually a matter *
! of using the binomial theorem. This is what we do below and to make  *
! the computation efficient the expression (2.4) is written as a       *
! succint double sum.                                                  *
!----------------------------------------------------------------------*
      Subroutine fFactor(loneX,ltwoX,lsumX,loneY,ltwoY,lsumY,loneZ      &
     &                 ,ltwoZ,lsumZ,PAxyz,PBxyz,FactorX,FactorY,FactorZ)
      Implicit Real*8 (a-h,o-z)

      Parameter(MaxAngqNr=6,MaxAr=MaxAngqNr*(MaxAngqNr+1)/2)
      Dimension FactorX(2*MaxAngqNr+1),FactorY(2*MaxAngqNr+1)
      Dimension FactorZ(2*MaxAngqNr+1)
      Dimension PAxyz(3),PBxyz(3)

      Do 105, ia=0,lsumX !We use unrolled loops with regard to x,y and z
        fff2=0          !therefore, here we start with the x-factors.
        iLowB=max(0,ia-ltwoX)  !These lower and upper bounds have to do
        iUpB=min(ia,loneX)   !with the allowed numbers in the binomial
        Do 106, i=iLowB,iUpB  !coefficients.
          fff1=NoverP_Q(loneX,i)*NoverP_Q(ltwoX,ia-i)
          If(i.ne.0) then  !This is needed for some compilers (NAG_64)
            PAraise=PAxyz(1)**i
          Else
            PAraise=1.0d0
          Endif
          If(ia-i.ne.0) then
            PBraise=PBxyz(1)**(ia-i)
          Else
            PBraise=1.0d0
          Endif
          fff2=fff2+fff1*PAraise*PBraise
106     Continue
        FactorX(lsumX-ia+1)=fff2
105   Continue
      Do 115, ia=0,lsumY  !y-factors.
        fff2=0
        iLowB=max(0,ia-ltwoY)
        iUpB=min(ia,loneY)
        Do 116, i=iLowB,iUpB
          fff1=NoverP_Q(loneY,i)*NoverP_Q(ltwoY,ia-i)
          If(i.ne.0) then  !This is needed for some compilers (NAG_64)
            PAraise=PAxyz(2)**i
          Else
            PAraise=1.0d0
          Endif
          If(ia-i.ne.0) then
            PBraise=PBxyz(2)**(ia-i)
          Else
            PBraise=1.0
          Endif
          fff2=fff2+fff1*PAraise*PBraise
116     Continue
        FactorY(lsumY-ia+1)=fff2
115   Continue
      Do 125, ia=0,lsumZ  !z-factorz.
        fff2=0
        iLowB=max(0,ia-ltwoZ)
        iUpB=min(ia,loneZ)
        Do 126, i=iLowB,iUpB
          fff1=NoverP_Q(loneZ,i)*NoverP_Q(ltwoZ,ia-i)
          If(i.ne.0) then  !This is needed for some compilers (NAG_64)
            PAraise=PAxyz(3)**i
          Else
            PAraise=1.0d0
          Endif
          If(ia-i.ne.0) then
            PBraise=PBxyz(3)**(ia-i)
          Else
            PBraise=1.0d0
          Endif
          fff2=fff2+fff1*PAraise*PBraise
126     Continue
        Factorz(lsumZ-ia+1)=fff2
125   Continue
      Return
      End
