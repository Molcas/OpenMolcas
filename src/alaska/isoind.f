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
      subroutine isoind(iSD,nSD,ish,nIrrep,IrrCmp,MxUnq,IndS,MxShll,
     &                  iAOtSO,MxAO)
      Integer iSD(0:nSD,1024), iTwoj(0:7), IrrCmp(MxUnq),
     &        IndS(MxShll), iAOtSO(MxAO,0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*
*...  generates SAO labels for a given shell
*
      icmp=iSD(2,ish)
      ibas=iSD(3,ish)
      iAO =iSD(7,ish)
      iSheli=iSD(11,ish)
*
      in2 = 0
      Do 10 j = 0, nIrrep-1
      Do 11 i1 = 1, iCmp
         If (iAnd(IrrCmp(IndS(iSheli)+i1),
     &             iTwoj(j)).ne.0) Then

         iSO = iAOtSO(iAO+i1,j)
         Do 12 iAOi = 0, iBas-1
            iSOi = iSO + iAOi
            in2 = in2 + 1
            write(6,*) 'Shell=',ish,'  in2=',in2,' Ang. component=',
     >                 i1,'  Contraction',iAOi+1,'  SO=',isoi,j+1
 12      Continue
         End If
 11   Continue
 10   Continue
      return
      end
*
