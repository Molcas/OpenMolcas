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
        subroutine JoinLvec (wrk,wrksize,
     c             PossV1,PossV2,NaGrpR,LunAux)
c
c       This routine do:
c        Join L0-L2 files (Global, dimensioned as nc)
c        from local L0-L2 files (dimensioned as ncLoc)
c        N.B. This file have sense only for parallel run
c
c
c       Structure of Cholesky vector files
c
c       L0(m,IJ)    L0vctr  I>=J
c
c       L1(m,I ,A') L1vcxx xx - Group of A'
c@@        kokot som, ze som to takto urobil, prerobit na L(m,a,i) to treba
c
c       L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
c                          yy - Group of B'
c
c       Memory requirements:
c        real:
c       V1   - max {ov'm; v'v'm;  oom}
c       V2   - max {v'v'm, ov'm; oom}
c        actual:
c       reord routine requirements are used (DistMemReord)
c
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
c
c       help variables
        character*6 LunName
        integer aGrp,bGrp
        integer dim1,dima,dimb,dimab
        integer i,ncLoc,ncOff
c
c
c7.0        calculate ncOffset for given node
        ncLoc=NChLoc(myRank)
        ncOff=0
        if (myRank.gt.0) then
          do i=0,myRank-1
          ncOff=ncOff+NChLoc(i)
          end do
        end if
c
c
c
c7.1        Make global L0
        LunName=L0Name
c
c7.1.1  Read Local L0 = V2(ml,ij) on proper place in V2
        dim1=ncLoc*no*(no+1)/2
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c7.1.2        Vanish V1(m,ij)
        dim1=nc*no*(no+1)/2
        call mv0zero (dim1,dim1,wrk(PossV1))
c
c7.1.3        Insert V2(ml,ij) -> V1(m,ij)
        dim1=no*(no+1)/2
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
c
c##        Synchronizacny bod
c7.1.4        Allreduce V1
        dim1=nc*no*(no+1)/2
        call gadgop (wrk(PossV1),dim1,'+')
c
c7.1.5        Save L0 (Global)
        dim1=nc*no*(no+1)/2
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c
c
c7.2    make global L1 files
c
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
        LunName=L1Name(aGrp)
c
c7.2.1  read L1 = V2(ml,i,a') back into file
        dim1=no*dima*ncLoc
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c7.2.2        Vanish V1(m,i,a')
        dim1=nc*no*dima
        call mv0zero (dim1,dim1,wrk(PossV1))
c
c7.2.3        Insert V2(ml,i,a') -> V1(m,i,a')
        dim1=no*dima
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
c
c##        Synchronizacny bod
c7.2.4        Allreduce V1
        dim1=nc*no*dima
        call gadgop (wrk(PossV1),dim1,'+')
c
c7.2.5        Save V1(m,i,a') (Global)
        dim1=nc*no*dima
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
        end do
c
c
c
c2.3    Make glogal L2(m,a'b') files
c
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
c
c
c7.3.1  read L2 = V2(ml,a'b') back into file
        dim1=dimab*ncLoc
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c7.3.2        Vanish V1(m,a'b')
        dim1=nc*dimab
        call mv0zero (dim1,dim1,wrk(PossV1))
c
c7.3.3        Insert V2(ml,a',b') -> V1(m,a',b')
        dim1=dimab
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
c
c##        Synchronizacny bod
c7.3.4        Allreduce V1
        dim1=nc*dimab
        call gadgop (wrk(PossV1),dim1,'+')
c
c7.3.5        Save V1(ml,a',b') (Global)
        dim1=nc*dimab
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
        end do
        end do
c
#else
c Avoid unused argument warnings
        if (.false.) then
          call Unused_real_array(wrk)
          call Unused_integer(PossV1)
          call Unused_integer(PossV2)
          call Unused_integer(NaGrpR)
          call Unused_integer(LunAux)
        end if
#endif
c
        return
        end
