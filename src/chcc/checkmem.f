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
        subroutine checkMem(NvGrp, NvSGrp, NchBlk,
     & Jal1, Jal2, wrksize, maxdim)
c
c Check memory consumption according to the specified
c orbital segmentation
c
        implicit none
#include "chcc1.fh"
c
        integer NvGrp,NvSGrp,NchBlk
        integer wrksize
c
        integer maxdim,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
        integer PossV1,PossV2,PossV3,PossV4
        integer PossL11,PossL12
        integer PossH1,PossH2,PossH3,PossH4,PossH5
        integer PossM1,PossM2,PossM3,PossM4,PossM5
        integer PossK,PossQ
        integer PossT,PossMax
c        jalove
        integer Jal1,Jal2,Jal3
c
c2      Distribute memory (and redefine wrksize)
c
c*.1        Distribute Permanent arrays
c
        PossT=1
        call DistMemPerm (PossT)
        PossMax=PossT
c
c*.2    Distribute Work arrays for Reord
c
        call DefParReord (NvGrp,maxdim)
c        also
        call DefParo2v4 (NvGrp,NvGrp,NvSGrp,NvSGrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
        PossT=PossFree
        call DistMemReord (NvGrp,maxdim,mdSGrpa,NchBlk,
     c       PossV1,PossV2,PossV3,PossV4,PossM1,PossM2,
     c       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
c
c*.3    Distribute Work arrays for o3v3
c
        call DefParo3v3 (NvGrp,maxdim)
c
        PossT=PossFree
        call DistMemo3v3jk (NvGrp,maxdim,
     c       PossV1,PossV2,PossV3,PossV4,
     c       PossH1,PossH2,PossH3,PossH4,PossH5,
     c       PossK,PossQ,
     c       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
c
        PossT=PossFree
        call DistMemo3v3t2 (NvGrp,maxdim,
     c       PossV1,PossV2,PossV3,PossV4,
     c       PossH1,PossH2,PossH3,PossH4,
     c       PossK,PossQ,
     c       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
c
        PossT=PossFree
        call DistMemo3v3chol (NvGrp,maxdim,
     c       PossV1,PossV2,PossV3,PossV4,
     c       PossH1,PossH2,PossH3,PossH4,
     c       PossM1,PossM2,PossM3,PossM4,PossM5,
     c       PossK,PossQ,
     c       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
c
c*.4    Distribute Work arrays for o2v4
c
        call DefParo2v4 (NvGrp,NvGrp,NvSGrp,NvSGrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
        PossT=PossFree
c        Nazvy premennych tu davam ine, lebo ide len o vypocet narokov na Mem
        call DistMemo2v4 (NvGrp,NvGrp,NvSGrp,NvSGrp,
     c                 mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,
     c                 PossV1,PossV2,PossV3,PossV4,
     c                 PossL11,PossL12,
     c                 PossH1,PossH2,PossH3,PossH4,PossQ,
     c                 Jal1,Jal2,
     c                 PossM1,PossM2,PossM3,PossM4,PossM5,Jal3,
     c                 PossT,PossK)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
c
c*.5        Memory requirements for Summary step are upperbounded by o3v3
c        and also o2v4, this it is skipped
c
        wrksize=PossMax
c
        return
        end
