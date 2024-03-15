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
* Copyright (C) 2022, Jie J. Bao                                       *
************************************************************************
******************************************************************
* history:                                                       *
* Jie J. Bao, on Apr. 11, 2022, created this file.               *
******************************************************************

* This file contains subroutines relating to generalized 1-e
* density matrix (GD) called in CMSNewton, including
* CalcGD:       calculating GD with lucia.
* CalcDg:       calculating Dg matrix, namely sum_{vx}{GD^KL_vx * g_tuvx}
* RotGD:        GD^KL_tu = sum_{MN}{U^KM * U^LN * GD^MN_tu}
* TransposeMat: transform GD^KL_tu from leading with state indices
*               to leading with orbital indices in mode 1, and vice
*               versa in mode 2.



      Subroutine RotGD(GD,R,nGD,lRoots,NAC2)
      use CMS, only: RGD
      INTEGER nGD,lRoots,NAC2
      Real*8 GD(nGD),R(lRoots**2)

      INTEGER iNAC2,iLoc,lRoots2
*      Real*8 RGD(lRoots**2),RGDR(lRoots**2)

      lRoots2=lRoots**2


C      write(6,*) 'rotation matrix in RotGD'
C      CALL RecPrt(' ',' ',R,lRoots,lRoots)
C
C      write(6,*) 'GD matrix after rotation'
C      CALL RecPrt(' ',' ',GD,lRoots2,NAC2)


      DO iNAC2=1,NAC2
       iLoc=(iNAC2-1)*lRoots2+1
       CALL DGEMM_('T','N',lRoots,lRoots,lRoots,
     &               1.0d0,R       ,lRoots,GD(iLoc),lRoots,
     &               0.0d0,RGD     ,lRoots)
       CALL DGEMM_('N','N',lRoots,lRoots,lRoots,
     &               1.0d0,RGD     ,lRoots,R       ,lRoots,
     &               0.0d0,GD(iLoc),lRoots)
      END DO

C      write(6,*) 'GD matrix after rotation'
C      CALL RecPrt(' ',' ',GD,lRoots2,NAC2)

      RETURN
      End Subroutine


      Subroutine TransposeMat(Matout,Matin,nElem,nRow_in,nCol_in)
      INTEGER nElem,nRow_in,nCol_in,iRow,iCol,iOff1,iOff2
      Real*8 Matin(nElem),Matout(nElem)

      IF(nRow_in*nCol_in.ne.nElem) THEN
       write(6,*) 'Error in TransposeMat()'
       write(6,*) 'nRow_in*nCol_in != nElem'
      END IF

      DO iCol=1,nCol_in
       iOff1=(iCol-1)*nRow_in
       Do iRow=1,nRow_in
        iOff2=(iRow-1)*nCol_in
        Matout(iOff2+iCol)=Matin(iOff1+iRow)
       End Do
      END DO

      RETURN
      End Subroutine
************************************************************************

      Subroutine CalcDg(Dgorbit,GDorbit,Gtuvx,nGD,nTUVX,NAC,lRoots)
      INTEGER nGD,nTUVX,NAC,lRoots
      Real*8 Dgorbit(nGD),GDorbit(nGD),Gtuvx(nTUVX)

      INTEGER NAC2,lRoots2

      NAC2=NAC**2
      lRoots2=lRoots**2

      CALL DGEMM_('T','N',NAC2,lRoots2,NAC2,1.0d0,
     &              Gtuvx,NAC2,GDorbit,NAC2,0.0d0,
     &                         Dgorbit,NAC2)

      RETURN
      End Subroutine
************************************************************************

      Subroutine CalcGD(GD,nGD)
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_lucia, only: DStmp, Dtmp
      use Lucia_Interface, only: Lucia_Util
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      INTEGER nGD
      Real*8 GD(nGD)
      INTEGER CIDisk1,CIDisk2
      INTEGER p,q,ipq,iqp,NAC2,IOffNIJ1,IOffNIJ2
      Real*8, Allocatable:: SDtmp(:), TmpD(:)
      Real*8, Allocatable, Target:: VecL(:), VecR(:)

      NAC2=NAC**2
      Call mma_allocate(VecL,NConf,Label='VecL')
      Call mma_allocate(VecR,NConf,Label='VecR')
      Call mma_allocate(TmpD,NAC**2,Label='TmpD')
      Call mma_allocate(SDtmp,NAC**2,Label='SDtmp')
      SDtmp(:)=DStmp(:)
      TmpD(:)=DTmp(:)
      CIDisk1=IADR15(4)
      Do jRoot=1,lRoots
       Call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
       CIDisk2=IADR15(4)
       Do kRoot=1,jRoot-1
        Call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
        Call Lucia_Util('Densi',
     &                   CI_Vector=VecL(:),
     &                   RVec=VecR(:))
        IOffNIJ1=(lRoots*(jRoot-1)+kRoot-1)*NAC2
        IOffNIJ2=(lRoots*(kRoot-1)+jRoot-1)*NAC2
C        write(6,*)'GD matrix',jRoot,kRoot
C        CALL RecPrt(' ',' ',Dtmp,NAC,NAC)
        Call DCopy_(NAC2,Dtmp,1,GD(IOffNIJ1+1),1)
         dO q=1,NAC
          do p=1,NAC
          ipq=(q-1)*NAC+p
          iqp=(p-1)*NAC+q
          GD(IOffNIJ2+iqp)=Dtmp(ipq)
*          GDMat(NIJ2,q,p)=Dtmp(q+(p-1)*NAC)
          end do
         eND dO
       End Do
       kRoot=jRoot
       Call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
       Call Lucia_Util('Densi',
     &                 CI_Vector=VecL(:),
     &                 RVec=VecR(:))
       IOffNIJ1=(lRoots+1)*(jRoot-1)*NAC2
C       write(6,*)'GD matrix',jRoot,kRoot
C       CALL RecPrt(' ',' ',Dtmp,NAC,NAC)
       Call DCopy_(NAC2,Dtmp,1,GD(IOffNIJ1+1),1)
      End DO
      DStmp(:)=SDtmp(:)
      Dtmp(:)=TmpD(:)
      Call mma_deallocate(SDtmp)
      Call mma_deallocate(TmpD)
      Call mma_deallocate(VecL)
      Call mma_deallocate(VecR)
      END Subroutine
************************************************************************


