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
        subroutine distnodes
c
c        this routine distribute nodes to different parts
c
        use Para_Info, only: nProcs
        implicit none
#include "parallel.fh"
c
        integer i
        REAL*8 efftot
c
ctmp ta zatial takto
        if (nProcs.eq.1) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=0
        idaabb=0
        idabba=0
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.2) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=1
        idaabb=1
        idabba=1
c
cIII        def node for finale
        idfin=1
c
c
        else if (nProcs.eq.3) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=2
c
cIII        def node for finale
        idfin=1
c
c
        else if (nProcs.eq.4) then
c
cI        def nodes for sumoverab
        nprocab=4
        idab(1)=0
        idab(2)=1
        idab(3)=2
        idab(4)=3
        ideffab(1)=0.25
        ideffab(2)=0.25
        ideffab(3)=0.25
        ideffab(4)=0.25
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=3
c
cIII        def node for finale
        idfin=3
c
c
        else if (nProcs.eq.5) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=3
        idabba=4
c
cIII        def node for finale
        idfin=2
c
c
        else if (nProcs.eq.6) then
c
cI        def nodes for sumoverab
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
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=4
        idabba=5
c
cIII        def node for finale
        idfin=3
c
        else if (nProcs.eq.10) then
c
cI        def nodes for sumoverab
        nprocab=4
          idab(1)=0
          idab(2)=1
          idab(3)=2
          idab(4)=3
          ideffab(1)=1.0
          ideffab(2)=1.0
          ideffab(3)=1.0
          ideffab(4)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=4
        idbaab=5
        idbbaa=6
        idbbbb=7
        idaabb=8
        idabba=9
c
cIII        def node for finale
        idfin=5
c
        else
c
cI        def nodes for sumoverab
        nprocab=nProcs
        do i=1,nprocab
          idab(i)=i-1
          ideffab(i)=1.0
        end do
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=4
        idabba=5
c
cIII        def node for finale
        idfin=6

        end if
c
         return
c
ctmp         koniec tmp riesenia
c
c
c
        if (nProcs.eq.1) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=0
        idaabb=0
        idabba=0
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.2) then
c
cI        def nodes for sumoverab
        nprocab=2
        idab(1)=0
        idab(2)=1
        ideffab(1)=0.5
        ideffab(2)=0.5
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=1
        idaabb=1
        idabba=1
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.3) then
c
cI        def nodes for sumoverab
        nprocab=3
        idab(1)=0
        idab(2)=1
        idab(3)=2
        ideffab(1)=0.333
        ideffab(2)=0.333
        ideffab(3)=0.333
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=2
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.4) then
c
cI        def nodes for sumoverab
        nprocab=2
        idab(1)=2
        idab(2)=3
        ideffab(1)=1.0
        ideffab(2)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=1
        idaabb=1
        idabba=1
c
cIII        def node for finale
        idfin=0
c
c
        else
c
cI        def nodes for sumoverab
c
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
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=0
        idbbaa=0
        idbbbb=2
        idaabb=3
        idabba=3
c
cIII        def node for finale
        idfin=1
c
        end if
c
c        renormalize ideffab
c
        efftot=0.0d0
        do i=1,nprocab
        efftot=efftot+ideffab(i)
        end do
c
        do i=1,nprocab
        ideffab(i)=ideffab(i)/efftot
        end do

c
        return
        end
