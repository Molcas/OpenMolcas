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
! MP2 density correction.
      Subroutine DCorrCorr(Dens,DenCorr,Trace_Diff,iOrb,iOcc)
      Implicit Real*8 (a-h,o-z)
      Dimension Dens(*),DenCorr(*)
      Trace_HF=dble(iOcc*2)
      kaunt=0
      T=Trace_HF/(Trace_HF-Trace_Diff)
      Do 183, i=1,iOrb
        Do 184, j=1,i
          kaunt=kaunt+1
          Dens(kaunt)=T*(Dens(kaunt)-DenCorr(kaunt))
184     Continue
183   Continue
!      Trace=0.0d0
!      kaunt=0
!      Do 181, i=1,iOrb
!        Do 182, j=1,i
!          kaunt=kaunt+1
!          If(i.eq.j)Trace=Trace+Dens(kaunt)
!182     Continue
!181   Continue
!      call triprt('KKK',' ',Dens,iorb)
!      write(6,*)'QQQ:',Trace
      Return
      End
