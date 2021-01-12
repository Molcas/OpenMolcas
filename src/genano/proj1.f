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
* The procedure is to read the projection orbitals, generate a         *
* complimentary space, and orthonormalize it. Then the density matrix  *
* is transformed into this representation, the corresponding matrix    *
* are zeroed out, and transformed back to AO representaion.            *
*                                                                      *
*======================================================================*
*                                                                      *
* NOTE: The transformations and orthonormalizations are done as        *
*       (On^4) processes and not as O(n^3) processed due to a lazy     *
*       programmer. I will rewrite this when I get the time.           *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine Proj1
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
      Dimension pOrb(MxS,MxS)
      Dimension DAO(MxS,MxS)
      Dimension DMO(MxS,MxS)
      Dimension Sc(MxS,MxS)
      Dimension Ovlp(MxS,MxS)
      call molcas_open(17,'PROJ')
c      Open(unit=17,file='PROJ',form='FORMATTED')
      iBlk=0
*--- Loop over l quantum number ---*
      Do 100 iLqn=0,MxLqn
         Read(17,*,end=999,err=999) nB,nO
         If(nB*nO.le.0) Go To 100
         If(nB.ne.nPrim(iLqn)) Then
            Write(6,*) 'Project: inconsistency in number of functions'
            Call Abend
         End If
         Do 105 iB=1,nB
            Read(17,*) (pOrb(iB,iO),iO=1,nO)
105      Continue
*--- Generate complimentary space ---*
         Do 110 iB=1,nB
         Do 111 iO=nO+1,nB
            pOrb(iB,iO)=Cos(3.14d0*iO*iB/nB)
111      Continue
110      Continue
*         Write(*,'(a,i5)') ' Projection orbitals',iLqn
*         Do 120 iB=1,nB
*            Write(*,'(1x,6f12.6)') (pOrb(iB,iO)),iO=1,nO)
*120      Continue
*--- Loop over m quantum number ---*
         Do 200 iShell=-iLqn,iLqn
            iBlk=iBlk+1
*           Write(*,'(a,2i5)') ' Density matrix block',iLqn,iShell
*           Call Triprt(' ','(6F12.6)',tDsym(iSymbk(iBlk)),nPrim(iLqn))
*           Write(*,'(a,2i5)') ' Overlap matrix block',iLqn,iShell
*           Call Triprt(' ','(6F12.6)',Ssym(iSymbk(iBlk)),nPrim(iLqn))
*--- Copy and square D and S ---*
            Do 210 i=1,nB
            Do 211 j=1,nB
               indT=min(i,j)+max(i,j)*(max(i,j)-1)/2
               DAO(i,j)=tDsym(iSymbk(iBlk)-1+indT)
               Ovlp(i,j)=Ssym(iSymbk(iBlk)-1+indT)
211         Continue
210         Continue
*--- Orthonormalize orbitals ---*
            Do 220 iO=1,nB
               s=0.0d0
               Do 221 iB=1,nB
               Do 2210 jB=1,nB
                  s=s+pOrb(iB,iO)*ovlp(iB,jB)*pOrb(jB,iO)
2210           Continue
221            Continue
               s=1.0d0/sqrt(s)
               Do 222 iB=1,nB
                  pOrb(iB,iO)=s*pOrb(iB,iO)
222            Continue
               Do 223 jO=iO+1,nB
                  s=0.0d0
                  Do 224 iB=1,nB
                  Do 2240 jB=1,nB
                     s=s+pOrb(iB,iO)*ovlp(iB,jB)*pOrb(jB,jO)
2240              Continue
224               Continue
                  Do 225 iB=1,nB
                     pOrb(iB,jO)=pOrb(iB,jO)-s*pOrb(iB,iO)
225               Continue
223            Continue
220         Continue
*--- Transform to MO basis ---*
            Do 230 iB=1,nB
            Do 231 iO=1,nB
               s=0.0d0
               Do 232 jB=1,nB
                  s=s+ovlp(iB,jB)*pOrb(jB,iO)
232            Continue
               Sc(iB,iO)=s
231         Continue
230         Continue
            Do 240 iO=1,nB
            Do 241 jO=1,nB
               s=0.0d0
               Do 242 iB=1,nB
               Do 243 jB=1,nB
                  s=s+Sc(iB,iO)*DAO(iB,jB)*Sc(jB,jO)
243            Continue
242            Continue
               DMO(iO,jO)=s
241         Continue
240         Continue
*--- Project ---*
            Do 250 iO=1,nO
            Do 251 jO=1,nB
               DMO(iO,jO)=0.0d0
               DMO(jO,iO)=0.0d0
251         Continue
250         Continue
*--- Transform back to AO ---*
            ind=0
            Do 260 iB=1,nB
            Do 261 jB=1,iB
               ind=ind+1
               s=0.0d0
               Do 262 iO=1,nB
               Do 263 jO=1,nB
                  s=s+pOrb(iB,iO)*DMO(iO,jO)*pOrb(jB,jO)
263            Continue
262            Continue
               tDsym(iSymbk(iBlk)-1+ind)=s
261         Continue
260         Continue
*--- ---*
*           Write(*,'(a,2i5)') ' Density matrix block',iLqn,iShell
*           Call Triprt(' ','(6F12.6)',tDsym(iSymbk(iBlk)),nPrim(iLqn))
200      Continue
100   Continue
*--- Finished ---*
999   Continue
      Close(unit=17)
      Return
      End
