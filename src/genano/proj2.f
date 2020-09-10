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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
* This routine projects out specified vectors from the density matrix. *
* The procedure is to read the projection orbitals and perform a       *
* sequence of rank-1 updated to zero out the contributions <c|D|c>.    *
* This procedure does not necessarily produce a complete projection    *
* of the density matrix, and routine 'proj1' is preferable.            *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine Proj2
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
      Dimension ProjOrb(MxS*MxS)
      Dimension TmpDens(MxS*MxS)
      Dimension TmpOvlp(MxS*MxS)
      call molcas_open(17,'PROJ')
c      Open(unit=17,file='PROJ',form='FORMATTED')
      iBlk=0
*--- Loop over l quantum number ---*
      Do 100 iLqn=0,MxLqn
*--- Read projection orbitals ---*
         Read(17,*,end=999,err=999) nB,nO
         If(nB*nO.le.0) Go To 100
         If(nB.ne.nPrim(iLqn)) Then
            Write(6,*) 'Project: inconsistency in number of functions'
            Call Abend
         End If
         Do 105 iB=1,nB
            Read(17,*) (ProjOrb(iB+nB*(iO-1)),iO=1,nO)
105      Continue
*         Write(*,'(a,i5)') ' Projection orbitals',iLqn
*         Do 101 iB=1,nB
*            Write(*,'(1x,6f12.6)') (ProjOrb(iB+nB*(iO-1)),iO=1,nO)
*101      Continue
*--- Loop over m quantum number ---*
         Do 110 iShell=-iLqn,iLqn
            iBlk=iBlk+1
*           Write(*,'(a,2i5)') ' Density matrix block',iLqn,iShell
*           Call Triprt(' ','(6F12.6)',tDsym(iSymbk(iBlk)),nPrim(iLqn))
*           Write(*,'(a,2i5)') ' Overlap matrix block',iLqn,iShell
*           Call Triprt(' ','(6F12.6)',Ssym(iSymbk(iBlk)),nPrim(iLqn))
*--- Copy S and D into square form ---*
            Do 111 i=1,nB
            Do 1110 j=1,nB
               indT=min(i,j)+max(i,j)*(max(i,j)-1)/2
               indS=i+nB*(j-1)
               TmpDens(indS)=tDsym(iSymbk(iBlk)-1+indT)
               TmpOvlp(indS)=Ssym(iSymbk(iBlk)-1+indT)
*              Write(*,'(1x,4i5,2f12.6)') i,j,IndT,IndS,
*    &            TmpDens(indS),TmpOvlp(indS)
1110        Continue
111         Continue
*--- Project ---*
            Do 112 iO=1,nO
               eval=0.0d0
               Do 113 iB=1,nB
               Do 1130 jB=1,nB
               Do 1131 kB=1,nB
               Do 1132 lB=1,nB
                  ij=iB+nB*(jB-1)
                  jk=jB+nB*(kB-1)
                  kl=kB+nB*(lB-1)
                  t=ProjOrb(iB+nB*(iO-1))*ProjOrb(lB+nB*(iO-1))
                  t=t*TmpOvlp(ij)
                  t=t*TmpDens(jk)
                  t=t*TmpOvlp(kl)
                  eval=eval+t
1132           Continue
1131           Continue
1130           Continue
113            Continue
               cnorm=0.0d0
               Do 114 iB=1,nB
               Do 1140 jB=1,nB
                  ij=iB+nB*(jB-1)
                  t=ProjOrb(iB+nB*(iO-1))*ProjOrb(jB+nB*(iO-1))
                  t=t*TmpOvlp(ij)
                  cnorm=cnorm+t
1140           Continue
114            Continue
*              Write(*,'(a,f12.6)') ' Eigenvalue ',eval
*              Write(*,'(a,f12.6)') ' Norm       ',cnorm
               eval=eval/cnorm/cnorm
               Do 115 iB=1,nB
               Do 1150 jB=1,nB
                  t=ProjOrb(iB+nB*(iO-1))*ProjOrb(jB+nB*(iO-1))
                  TmpDens(iB+nB*(jB-1))=TmpDens(iB+nB*(jB-1))-eval*t
1150           Continue
115            Continue
112         Continue
*--- Copy back D ---*
            Do 118 i=1,nB
            Do 1180 j=1,nB
               indT=min(i,j)+max(i,j)*(max(i,j)-1)/2
               indS=i+nB*(j-1)
               tDsym(iSymbk(iBlk)-1+indT)=TmpDens(indS)
1180        Continue
118         Continue
*           Write(*,'(a,2i5)') ' Density matrix block',iLqn,iShell
*           Call Triprt(' ','(6F12.6)',tDsym(iSymbk(iBlk)),nPrim(iLqn))
110      Continue
100   Continue
*--- ---*
999   Continue
      Close(unit=17)
      Return
      End
