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
* Copyright (C) 1997, Anders Bernhardsson                              *
************************************************************************
      Subroutine Niclas(H,coor,LUT)
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iChTbl
* eaw 970909
      Implicit Real*8(a-h,o-z)
#include "stdalloc.fh"
#include "SysDef.fh"
      Real*8 H(*)
      Character*40 Label
      Integer nDeg(200),ldisp(0:7)
      Integer inddsp(100,0:7)
      Logical, External :: TF
      Real*8 Coor(*)
      Real*8 Dummy(1)
      Real*8, Allocatable:: Htmp(:), Tmp(:)
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      irec(i,j)=nd*(j-1)+i-1
*
      idsp=0
      Call iCOPY(nirrep,[0],0,ldisp,1)
      Do iIrrep=0,nIrrep-1
      mdc=0
       Do iCnttp = 1, nCnttp
        nCnti = dbsc(iCnttp)%nCntr
        Do iCnt = 1, nCnti
         mdc=mdc+1
         IndDsp(mdc,iIrrep)=idsp
         Do iCar = 0, 2
          iComp = 2**iCar
          If (TF(mdc,iIrrep,iComp)) Then
             idsp=idsp+1
             ldisp(iirrep)=ldisp(iirrep)+1
             ndeg(idsp)=nIrrep/dc(mdc)%nStab
          End If
         End Do
        End Do
       End Do
      End Do
*
************************************************************************
*
*    Steady
*
*    Make the symmetrized Hessian correct for degenerated geometries
*
************************************************************************
*
      nd=0
      Do i=0,nIrrep-1
       nD=ldisp(i)+nd
      End Do
      Call mma_allocate(TMP,nd**2,Label='Tmp')
      Call mma_allocate(HTMP,nd**2,Label='Htmp')
      Htmp(:)=0.0d0
      ii=0
      iii=0
      Do iS=1,Nirrep
       Do i = 1, ldisp(iS-1)
        Do j=1,i
            Tmp(itri(iii+i,iii+j)) =
     &                 sqrt(DBLE(nDeg(i+iii)*nDeg(j+iii)))*
     &                  H(ii+itri(i,j))
*          Write(*,*) H(ii+itri(i,j)),Tmo(itri(iii+i,iii+j))
        End Do
       End Do
       ii=ii+ldisp(is-1)*(ldisp(is-1)+1)/2
       iii=iii+ldisp(is-1)
      End Do
*
********************************************************************************
*
*   Go
*
********************************************************************************
*
      Call FCOOR(LUT,Coor)
      mdc=0
      iPERT=0
      Do iCnttp = 1, nCnttp
       nCnti = dbsc(iCnttp)%nCntr
       Do iCnt = 1, nCnti
        mdc=mdc+1
*
        nCenti=nIrrep/dc(mdc)%nStab
*
      ndc=0
      jPERT=0
      Do jCnttp = 1, nCnttp
       nCntj = dbsc(jCnttp)%nCntr
       Do jCnt = 1, nCntj
        ndc=ndc+1

        nCentj=nIrrep/dc(ndc)%nStab
        Do iIrrep=0,nIrrep-1
         iDsp = IndDsp(mdc,iIrrep)
         Do iCar = 0, 2
          iComp = 2**iCar
          If (TF(mdc,iIrrep,iComp)) Then
            idsp=idsp+1
            jDsp = IndDsp(ndc,iIrrep)
            Do jCar = 0, 2
             jComp = 2**jCar
             If (TF(ndc,iIrrep,jComp)) Then
              jdsp=jdsp+1
              HE=Tmp(itri(idsp,jdsp))
              Do iCo=0,Ncenti-1
               Do jCo=0,Ncentj-1
                i=iPert+ico*3+icar+1
                j=jPert+jco*3+jcar+1
                kop_m=dc(mdc)%iCoSet(iCo,0)
                nop_m=nropr(kop_m)
                kop_n=dc(ndc)%iCoSet(jCo,0)
                nop_n=nropr(kop_n)
                riPh=DBLE(iPrmt(nop_m,icomp)*iChTbl(iIrrep,nop_m))
     &           /sqrt(DBLE(nCENTI))
                rjPh=DBLE(iPrmt(nop_n,jcomp)*ichtbl(iirrep,nop_n))
     &          /sqrt(DBLE(nCENTJ))
                Htmp(1+irec(i,j))=Htmp(1+irec(i,j))+riph*rjph*HE
               End Do ! jco
              End Do ! ico
            End If
          End Do ! jcar
            End If
          End do ! icar
         End Do ! irrep
        jPert=jpert+ncentj*3
       End Do ! jcnt
      End Do ! jcnttp
        iPert=ipert+ncenti*3
       End Do ! icnt
      End Do ! icnttp
*
      Label='Unsymmetrized Hessian'
      WRITE(LUT,'(A)') Label
      Write(LUT,'(A)') '*BEGIN HESSIAN'
      Write(LUT,'(A,I5)') '*Number of pert. ',nd
      Call WRH(LUT,1,[nd],[nd],Htmp,Dummy,0,Label)
      Write(LUT,'(A)') '*END HESSIAN'
*
      Call Put_dArray('FC-Matrix',Htmp,nd**2)
*
      Call mma_deallocate(HTMP)
      Call mma_deallocate(TMP)
*
      Return
      End
