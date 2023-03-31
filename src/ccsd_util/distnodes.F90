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
        subroutine distnodes
!
!        this routine distribute nodes to different parts
!
        use Para_Info, only: nProcs
        implicit none
#include "parallel.fh"
!
        integer i
        REAL*8 efftot
!
!tmp ta zatial takto
        if (nProcs.eq.1) then
!
!I        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=0
        idaabb=0
        idabba=0
!
!III        def node for finale
        idfin=0
!
!
        else if (nProcs.eq.2) then
!
!I        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=1
        idaabb=1
        idabba=1
!
!III        def node for finale
        idfin=1
!
!
        else if (nProcs.eq.3) then
!
!I        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=2
!
!III        def node for finale
        idfin=1
!
!
        else if (nProcs.eq.4) then
!
!I        def nodes for sumoverab
        nprocab=4
        idab(1)=0
        idab(2)=1
        idab(3)=2
        idab(4)=3
        ideffab(1)=0.25
        ideffab(2)=0.25
        ideffab(3)=0.25
        ideffab(4)=0.25
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=3
!
!III        def node for finale
        idfin=3
!
!
        else if (nProcs.eq.5) then
!
!I        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=3
        idabba=4
!
!III        def node for finale
        idfin=2
!
!
        else if (nProcs.eq.6) then
!
!I        def nodes for sumoverab
        nprocab=6
        idab(1)=0
        idab(2)=1
        idab(3)=2
        idab(4)=3
        idab(5)=4
        idab(6)=5
        ideffab(1)=1.0
        ideffab(2)=1.0
        ideffab(3)=1.0
        ideffab(4)=1.0
        ideffab(5)=1.0
        ideffab(6)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=4
        idabba=5
!
!III        def node for finale
        idfin=3
!
        else if (nProcs.eq.10) then
!
!I        def nodes for sumoverab
        nprocab=4
          idab(1)=0
          idab(2)=1
          idab(3)=2
          idab(4)=3
          ideffab(1)=1.0
          ideffab(2)=1.0
          ideffab(3)=1.0
          ideffab(4)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=4
        idbaab=5
        idbbaa=6
        idbbbb=7
        idaabb=8
        idabba=9
!
!III        def node for finale
        idfin=5
!
        else
!
!I        def nodes for sumoverab
        nprocab=nProcs
        do i=1,nprocab
          idab(i)=i-1
          ideffab(i)=1.0
        end do
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=4
        idabba=5
!
!III        def node for finale
        idfin=6

        end if
!
         return
!
!tmp         koniec tmp riesenia
!
!
!
        if (nProcs.eq.1) then
!
!I        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=0
        idaabb=0
        idabba=0
!
!III        def node for finale
        idfin=0
!
!
        else if (nProcs.eq.2) then
!
!I        def nodes for sumoverab
        nprocab=2
        idab(1)=0
        idab(2)=1
        ideffab(1)=0.5
        ideffab(2)=0.5
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=1
        idaabb=1
        idabba=1
!
!III        def node for finale
        idfin=0
!
!
        else if (nProcs.eq.3) then
!
!I        def nodes for sumoverab
        nprocab=3
        idab(1)=0
        idab(2)=1
        idab(3)=2
        ideffab(1)=0.333
        ideffab(2)=0.333
        ideffab(3)=0.333
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=2
!
!III        def node for finale
        idfin=0
!
!
        else if (nProcs.eq.4) then
!
!I        def nodes for sumoverab
        nprocab=2
        idab(1)=2
        idab(2)=3
        ideffab(1)=1.0
        ideffab(2)=1.0
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=1
        idaabb=1
        idabba=1
!
!III        def node for finale
        idfin=0
!
!
        else
!
!I        def nodes for sumoverab
!
        nprocab=nProcs-2
        idab(1)=1
        idab(2)=2
        idab(3)=4
        idab(4)=5
        idab(5)=6
        idab(6)=7
        do i=1,nprocab
        ideffab(i)=1.0d0/nprocab
        end do
        ideffab(1)=ideffab(1)/2
        ideffab(2)=ideffab(2)/2
!
!II        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=0
        idbbaa=0
        idbbbb=2
        idaabb=3
        idabba=3
!
!III        def node for finale
        idfin=1
!
        end if
!
!        renormalize ideffab
!
        efftot=0.0d0
        do i=1,nprocab
        efftot=efftot+ideffab(i)
        end do
!
        do i=1,nprocab
        ideffab(i)=ideffab(i)/efftot
        end do

!
        return
        end
