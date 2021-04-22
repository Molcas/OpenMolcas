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
        subroutine GetChVHlp3 (L2,Tmp,cGrp,deGrp,LunAux)
c
c       this routine do:
c       read L2(m,c'de') into L2 from disk file
c       structure of files: for each c',de' one file with
c       name L2Name(cGrp,deGrp)
c       @ citanie zatial odflaknute
c
c       parameter description:
c       L2     - Array for L2 (O)
c       Tmp    - Temporary array of L2 size, used for mapping, if needed
c       xGrp   - Groups of c, delta (I)
c       LunAux - lun for auxiliary reading (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
c
        real*8 L2(1)
        real*8 Tmp(1)
        integer cGrp,deGrp,LunAux
        character*6 LunName
c
c       help variables
        integer length
c
c        nacitanie (+expanzia, ak treba)
        if (cGrp.gt.deGrp) then
          LunName=L2Name(cGrp,deGrp)
          length=nc*DimGrpa(cGrp)*DimGrpbe(deGrp)
          call GetX (L2(1),length,LunAux,LunName,1,1)
        else if (cGrp.eq.deGrp) then
          LunName=L2Name(cGrp,deGrp)
          length=nc*DimGrpa(cGrp)*(DimGrpbe(deGrp)+1)/2
          call GetX (Tmp(1),length,LunAux,LunName,1,1)
          length=DimGrpa(cGrp)*(DimGrpbe(deGrp)+1)/2
          call Exp1 (Tmp(1),L2(1),nc,length,DimGrpa(cGrp))
        else
          LunName=L2Name(deGrp,cGrp)
          length=nc*DimGrpa(cGrp)*DimGrpbe(deGrp)
          call GetX (Tmp(1),length,LunAux,LunName,1,1)
          call Map3_132 (Tmp(1),L2(1),nc,DimGrpbe(deGrp),DimGrpa(cGrp))
        end if
c
        return
        end
