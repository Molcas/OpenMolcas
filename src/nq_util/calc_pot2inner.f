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
* Copyright (C) 2021, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Dec. 25, 2021, created this file.               *
* ****************************************************************
      Subroutine Calc_Pot2_Inner(Pot2,mGrid,MOP,MOU,MOV,MOX,lSum)
#include "nq_info.fh"
#include "stdalloc.fh"
******Input
      INTEGER mGrid
      Logical lSum
******Note: when lSum is .true., calculate P(U'VX+UV'X+UVX'),
******      otherwise calculate PUVX
      Real*8,DIMENSION(mGrid*nOrbt)::MOP,MOU,MOV,MOX
******Output
      Real*8,DIMENSION(nPot2)::Pot2

******Intermediate
      Real*8,DIMENSION(:),Allocatable::MOUVX,MOVX1,MOVX2
      INTEGER iGrid,iOff0,iOff1,iOff2,iOff3,iStack,
     &        nnUVX,iVX,
     &        pIrrep,uIrrep,vIrrep,xIrrep,xMax,puIrrep,
     &        u,v,x,vorb,xorb,ioffu,nporb


      CALL mma_allocate(MOVX1,nVXt*mGrid)
      IF(lSum) CALL mma_allocate(MOVX2,nVXt*mGrid)
      CALL mma_allocate(MOUVX,nUVXt*mGrid)

      DO vIrrep=0,mIrrep-1
       DO xIrrep=0,vIrrep
        Do v=1,nAsh(vIrrep)
         vorb=v+OffOrb(vIrrep)+nIsh(vIrrep)
         IOff1=(vorb-1)*mGrid
         IF(xIrrep.eq.vIrrep) THEN
          xMax=v
          iStack=v*(v-1)/2
         ELSE
          xMax=nAsh(xIrrep)
          iStack=(v-1)*xMax
         END IF
         Do x=1,xMax
          xorb=x+OffOrb(xIrrep)+nIsh(xIrrep)
          IOff2=(xorb-1)*mGrid
          IOff3=(OffVX(xIrrep,vIrrep)+iStack+x-1)*mGrid
      do iGrid=1,mGrid
       MOVX1(iGrid+IOff3)=MOV(iGrid+IOff1)*MOX(iGrid+IOff2)
      end do
         End Do
        End Do
       END DO
      END DO


      IF(lSum) THEN
       DO vIrrep=0,mIrrep-1
        DO xIrrep=0,vIrrep
         Do v=1,nAsh(vIrrep)
          vorb=v+OffOrb(vIrrep)+nIsh(vIrrep)
          IOff1=(vorb-1)*mGrid
          If(xIrrep.eq.vIrrep) Then
           xMax=v
           iStack=v*(v-1)/2
          Else
           xMax=nAsh(xIrrep)
           iStack=(v-1)*xMax
          End if
          Do x=1,xMax
           xorb=x+OffOrb(xIrrep)+nIsh(xIrrep)
           IOff2=(xorb-1)*mGrid
           IOff3=(OffVX(xIrrep,vIrrep)+iStack+x-1)*mGrid
      do iGrid=1,mGrid
       MOVX1(iGrid+IOff3)=MOVX1(iGrid+IOff3)+
     &                    MOX(iGrid+IOff1)*MOV(iGrid+IOff2)
       MOVX2(iGrid+IOff3)=MOU(iGrid+IOff1)*MOV(iGrid+IOff2)
      end do
          End Do
         End Do
        END DO
       END DO
      END IF



      DO uIrrep=0,mIrrep-1
       IOffU=OffOrb(uIrrep)+nIsh(uIrrep)
       DO vIrrep=0,mIrrep-1
        DO xIrrep=0,vIrrep
         Do iVX=1,nVX(xIrrep,vIrrep)
          IOff1=(Offvx(xIrrep,vIrrep)+iVX-1)*mGrid
          IOff0=OffUVX(xIrrep,vIrrep,uIrrep)+(iVX-1)*nAsh(uIrrep)
          Do u=1,nAsh(uIrrep)
           IOff2=(iOffU+u-1)*mGrid
           IOff3=(IOff0+u-1)*mGrid
      do iGrid=1,mGrid
       MOUVX(iGrid+IOff3)=MOU(iGrid+IOff2)*MOVX1(iGrid+IOff1)
      end do
          End Do
         End Do
        End Do
       END DO
      END DO

      IF(lSum) THEN
       DO uIrrep=0,mIrrep-1
        IOffU=OffOrb(uIrrep)+nIsh(uIrrep)
        DO vIrrep=0,mIrrep-1
         DO xIrrep=0,vIrrep
          Do iVX=1,nVX(xIrrep,vIrrep)
           IOff1=(Offvx(xIrrep,vIrrep)+iVX-1)*mGrid
           IOff0=OffUVX(xIrrep,vIrrep,uIrrep)+(iVX-1)*nAsh(uIrrep)
           Do u=1,nAsh(uIrrep)
            IOff2=(iOffU+u-1)*mGrid
            IOff3=(IOff0+u-1)*mGrid
      do iGrid=1,mGrid
       MOUVX(iGrid+IOff3)=MOUVX(iGrid+IOff3)+
     &                    MOX(iGrid+IOff2)*MOVX2(iGrid+IOff1)
      end do
           End Do
          End Do
         End Do
        END DO
       END DO
      END IF

******Use dgemm to calculate PUVX at this grid point
      DO pIrrep=0,mIrrep-1
       nporb=mOrb(pIrrep)
       IF(nporb.eq.0) CYCLE
       IF(nAsh(pIrrep).eq.0) CYCLE
       IOff1=OffOrb(pIrrep)*mGrid+1
       IOff2=OffPUVX(pIrrep)+1
       DO uIrrep=0,mIrrep-1
        puIrrep=IEOR(pIrrep,uIrrep)
        DO vIrrep=0,mIrrep-1
         xIrrep=IEOR(puIrrep,vIrrep)
         nnUVX=nUVX(xIrrep,vIrrep,uIrrep)
       IF((xIrrep.gt.vIrrep).or.(nnUVX.eq.0)) CYCLE
         IOff3=OffUVX(xIrrep,vIrrep,uIrrep)*mGrid+1
       CALL DGEMM_('T','N',npOrb,nnUVX,mGrid,
     &             1.0d0,MOP(iOff1),mGrid,MOUVX(IOff3),mGrid,
     &             1.0d0,Pot2(iOff2),npOrb)
         IOff2=IOff2+nnUVX*npOrb
        END DO
       END DO
      END DO



      CALL mma_deallocate(MOVX1)
      IF(lSum) CALL mma_deallocate(MOVX2)
      CALL mma_deallocate(MOUVX)

      RETURN
      End Subroutine
