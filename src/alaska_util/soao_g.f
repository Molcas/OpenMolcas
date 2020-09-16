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
      Subroutine SOAO_g(iSD4,nSD,nSO, MemPrm, MemMax,
     &                  iBsInc,jBsInc,kBsInc,lBsInc,
     &                  iPrInc,jPrInc,kPrInc,lPrInc,
     &                  ipMem1,ipMem2, Mem1,  Mem2,
     &                  iPrint,iFnc, MemPSO)
      use Basis_Info, only: Shells
      Implicit Real*8 (a-h,o-z)
*
      Integer iSD4(0:nSD,4), iAnga(4), iCmpa(4),
     &        iFnc(4), iAO(4)
*
      iPrimi   = Shells(iSD4( 0,1))%nExp
      jPrimj   = Shells(iSD4( 0,2))%nExp
      kPrimk   = Shells(iSD4( 0,3))%nExp
      lPriml   = Shells(iSD4( 0,4))%nExp
      iBasi    = Shells(iSD4( 0,1))%nBasis
      jBasj    = Shells(iSD4( 0,2))%nBasis
      kBask    = Shells(iSD4( 0,3))%nBasis
      lBasl    = Shells(iSD4( 0,4))%nBasis
      Do iQuad = 1, 4
         iAnga(iQuad)  = iSD4( 1,iQuad)
         iCmpa(iQuad)  = iSD4( 2,iQuad)
         iAO(iQuad)    = iSD4( 7,iQuad)
      End Do
*
                  Call PSOAO1(nSO,MemPrm, MemMax,
     &                        iAnga, iCmpa, iAO  ,iFnc,
     &                        iBasi,iBsInc, jBasj,jBsInc,
     &                        kBask,kBsInc, lBasl,lBsInc,
     &                        iPrimi,iPrInc,jPrimj,jPrInc,
     &                        kPrimk,kPrInc,lPriml,lPrInc,
     &                        ipMem1,ipMem2, Mem1,  Mem2, MemPSO)
#ifdef _DEBUG_
         Write (6,*) ' ************** Memory partioning **************'
         Write (6,*) ' ipMem1=',ipMem1
         Write (6,*) ' Mem1=',Mem1
         Write (6,*) ' Mem2=',Mem2
         Write (6,*) ' Mem3=',Mem3
         Write (6,*) ' MemPSO=',MemPSO
         Write (6,*) ' MemScr=',MemScr
         Write (6,*) ' Mend=',Mend

         Write (6,*) ' iBasi,iBsInc=',iBasi,iBsInc
         Write (6,*) ' jBasj,jBsInc=',jBasj,jBsInc
         Write (6,*) ' kBasi,kBsInc=',kBask,kBsInc
         Write (6,*) ' lBasl,lBsInc=',lBasl,lBsInc
         Write (6,*) ' iPrimi,iPrInc=',iPrimi,iPrInc
         Write (6,*) ' jPrimj,jPrInc=',jPrimj,jPrInc
         Write (6,*) ' kPrimk,kPrInc=',kPrimk,kPrInc
         Write (6,*) ' lPriml,lPrInc=',lPriml,lPrInc
         Write (6,*) ' ***********************************************'
#endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iPrint)
      End
