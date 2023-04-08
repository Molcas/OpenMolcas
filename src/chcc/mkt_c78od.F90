!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
        subroutine MkT_C78od (T2,Tp,Tm,                                 &
     &             dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no)
!
!        this routine do:
!       T2n(be',ga',u,v) <<-
!        C7                    T2+(be',ga',uv)
!        C8                    T2-(be',ga',uv)
!        for beSGrp.ne.gaSGrp
!
        implicit none
        integer dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no
        real*8 T2(1:dimbe,1:dimga,1:no,1:no)
        real*8 Tp(1:dimbepp,1:dimgapp,1:no*(no+1)/2)
        real*8 Tm(1:dimbepp,1:dimgapp,1:no*(no-1)/2)
!
!        help variables
        integer u,v,be,ga,uvp,uvm,uup,bep,gap
!
!        diagonal part - contributoion from T+ only
        do u=1,no
        uup=u*(u+1)/2
!
!          cycle over be",ga"
          gap=addgapp
          do ga=1,dimgapp
            gap=gap+1
            bep=addbepp
            do be=1,dimbepp
              bep=bep+1
              T2(bep,gap,u,u)=T2(bep,gap,u,u)+Tp(be,ga,uup)
            end do
          end do
!
        end do
!
!
!        off diagonal - both T+,T-
        uvm=0
        do u=2,no
        uvp=u*(u-1)/2
        do v=1,u-1
          uvm=uvm+1
          uvp=uvp+1
!
!          cycle over be",ga"
          gap=addgapp
          do ga=1,dimgapp
            gap=gap+1
            bep=addbepp
            do be=1,dimbepp
              bep=bep+1
              T2(bep,gap,u,v)=T2(bep,gap,u,v)                           &
     &                       +Tp(be,ga,uvp)+Tm(be,ga,uvm)
              T2(bep,gap,v,u)=T2(bep,gap,v,u)                           &
     &                       +Tp(be,ga,uvp)-Tm(be,ga,uvm)
            end do
          end do
!
!          cycle over be",ga"
          gap=addgapp
          do ga=1,dimgapp
            gap=gap+1
            bep=addbepp
            do be=1,dimbepp
              bep=bep+1
!              T2(bep,gap,v,u)=T2(bep,gap,v,u)
!    c                       +Tp(be,ga,uvp)-Tm(be,ga,uvm)
            end do
          end do
!
        end do
        end do
!
!
        return
        end
