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
        subroutine JoinLvec (wrk,wrksize,                               &
     &             PossV1,PossV2,NaGrpR,LunAux)
!
!       This routine do:
!        Join L0-L2 files (Global, dimensioned as nc)
!        from local L0-L2 files (dimensioned as ncLoc)
!        N.B. This file have sense only for parallel run
!
!
!       Structure of Cholesky vector files
!
!       L0(m,IJ)    L0vctr  I>=J
!
!       L1(m,I ,A') L1vcxx xx - Group of A'
!@@        kokot som, ze som to takto urobil, prerobit na L(m,a,i) to treba
!
!       L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
!                          yy - Group of B'
!
!       Memory requirements:
!        real:
!       V1   - max {ov'm; v'v'm;  oom}
!       V2   - max {v'v'm, ov'm; oom}
!        actual:
!       reord routine requirements are used (DistMemReord)
!
#ifdef _MOLCAS_MPP_
        use Para_Info, only: MyRank
#endif
        implicit none
        integer PossV1,PossV2,NaGrpR,LunAux
#include "chcc1.fh"
#include "chcc_reord.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
#include "wrk.fh"
#ifdef _MOLCAS_MPP_
#include "parcc.fh"
!
!       help variables
        character*6 LunName
        integer aGrp,bGrp
        integer dim1,dima,dimb,dimab
        integer i,ncLoc,ncOff
!
!
!7.0        calculate ncOffset for given node
        ncLoc=NChLoc(myRank)
        ncOff=0
        if (myRank.gt.0) then
          do i=0,myRank-1
          ncOff=ncOff+NChLoc(i)
          end do
        end if
!
!
!
!7.1        Make global L0
        LunName=L0Name
!
!7.1.1  Read Local L0 = V2(ml,ij) on proper place in V2
        dim1=ncLoc*no*(no+1)/2
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
!
!7.1.2        Vanish V1(m,ij)
        dim1=nc*no*(no+1)/2
        call mv0zero (dim1,dim1,wrk(PossV1))
!
!7.1.3        Insert V2(ml,ij) -> V1(m,ij)
        dim1=no*(no+1)/2
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
!
!##        Synchronizacny bod
!7.1.4        Allreduce V1
        dim1=nc*no*(no+1)/2
        call gadgop (wrk(PossV1),dim1,'+')
!
!7.1.5        Save L0 (Global)
        dim1=nc*no*(no+1)/2
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!
!
!7.2    make global L1 files
!
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
        LunName=L1Name(aGrp)
!
!7.2.1  read L1 = V2(ml,i,a') back into file
        dim1=no*dima*ncLoc
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
!
!7.2.2        Vanish V1(m,i,a')
        dim1=nc*no*dima
        call mv0zero (dim1,dim1,wrk(PossV1))
!
!7.2.3        Insert V2(ml,i,a') -> V1(m,i,a')
        dim1=no*dima
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
!
!##        Synchronizacny bod
!7.2.4        Allreduce V1
        dim1=nc*no*dima
        call gadgop (wrk(PossV1),dim1,'+')
!
!7.2.5        Save V1(m,i,a') (Global)
        dim1=nc*no*dima
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
        end do
!
!
!
!2.3    Make glogal L2(m,a'b') files
!
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
        do bGrp=1,aGrp
        dimb=DimGrpaR(bGrp)
        if(aGrp.eq.bGrp) then
          dimab=dima*(dima+1)/2
        else
          dimab=dima*dimb
        end if
        LunName=L2Name(aGrp,bGrp)
!
!
!7.3.1  read L2 = V2(ml,a'b') back into file
        dim1=dimab*ncLoc
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
!
!7.3.2        Vanish V1(m,a'b')
        dim1=nc*dimab
        call mv0zero (dim1,dim1,wrk(PossV1))
!
!7.3.3        Insert V2(ml,a',b') -> V1(m,a',b')
        dim1=dimab
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
!
!##        Synchronizacny bod
!7.3.4        Allreduce V1
        dim1=nc*dimab
        call gadgop (wrk(PossV1),dim1,'+')
!
!7.3.5        Save V1(ml,a',b') (Global)
        dim1=nc*dimab
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
        end do
        end do
!
#else
! Avoid unused argument warnings
        if (.false.) then
          call Unused_real_array(wrk)
          call Unused_integer(PossV1)
          call Unused_integer(PossV2)
          call Unused_integer(NaGrpR)
          call Unused_integer(LunAux)
        end if
#endif
!
        return
        end
