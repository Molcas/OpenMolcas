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
      Subroutine PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,
     &                 maxlab,labs,PrTot,iPL,iPoint)
************************************************************************
*                                                                      *
*     purpose: printing of tables for different tensor properties      *
*                                                                      *
*     Short           logical option for either Short (Short           *
*                     =.true., the total electronic contribution)      *
*                     long (Short=.false., orbital contributions)      *
*                     output                                           *
*     sig             the sign factor for the electronic contri-       *
*                     bution, positive for ndiamagnetic shielding      *
*     nIrrep          the number of irreducible representations        *
*     nBas            the number of functions in each representa-      *
*     (0:nIrrep-1)    tion                                             *
*     nTot            the total number of elements supplied for        *
*                     each component of the property tensor; equal     *
*                     either to 1 (total electronic and total nuc-     *
*                     lear contributions) or to the dimension of       *
*                     the basis set.                                   *
*     Occ(1:nTot)     Occupation numbers for all eigenvectors,         *
*                     a dummy for Short outputs                        *
*     ThrSV           threshold for Occupation numbers; if             *
*                     Occ(i).le.ThrSV the contribution will not        *
*                     be printed                                       *
*     PrEl(1:nTot,    matrix elements for all components 1,2,...,      *
*          1:maxlab)  maxlab, nTot entries for each component          *
*                                                                      *
*     PrNu(1:maxlab)  nuclear contributions for each component         *
*     maxlab          total number of cartesian components             *
*     labs(1:maxlab)  labels for each component                        *
*     TotEl(1:6)      auxiliary storage area                           *
*     PrTot(1:maxlab) Total value for each component                   *
*     iPL             Print level                                      *
*     iPoint          The number of the center                         *
*                                                                      *
* 2000 Dept. of Chem. Phys., Univ. of Lund, Sweden                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Character*4 lbb4
      Character*16 labs(1:maxlab),lab
      Character*132 outlab
      Character*5 lab5
      Logical Short
      Integer nBas(0:nIrrep-1)
      Real*8 Occ(1:nTot), PrEl(1:nTot,1:maxlab),PrNu(1:maxlab),
     &       TotEl(1:6), PrTot(1:maxlab)
      Character Format1*18, Format2*13, Format3*14, Format4*14
*
      Format1='(2i5,f14.8,6f16.8)'
      Format2='(1x,a,6f16.8)'
      Format3='(1x,a,6f16.8/)'
      Format4='(a5,6f16.8)'
      Value=Zero
      Do i = 1, maxlab
         Value=Max(Value,Abs(PrNu(i)))
         Do j = 1, nTot
            Value=Max(Value,Abs(PrEl(j,i)))
         End Do
      End Do
      Value=Max(Value,One)
      nLog10=Max(1,INT(Log10(Value)+One)+1)
      nDec=Min(14-nLog10,8)
      Write(Format1(17:17),'(I1)') nDec
      Write(Format2(12:12),'(I1)') nDec
      Write(Format3(12:12),'(I1)') nDec
      Write(Format4(10:10),'(I1)') nDec
*
      If (.Not.Short) Write (6,'(A,A,D9.2/)')
     &' orbital contributions printed for occupation numbers',
     &' gt.',ThrSV
*
      Do i = 1, maxlab, 6
*
         Do j=1,132
            outlab(j:j)=' '
         End Do
*
         If (Short) Then
            If (iPL.ge.3) Then
              If (maxlab.gt.1.or.labs(1).ne.'                ') Then
                Write(outlab,'(1x,a,6a16)') 'Component              ',
     &                                   (labs(j),j=i,min(i+5,maxlab))
                Write(6,'(a)') outlab(:mylen(outlab))
              End If
            End If
            jcount=0
            Do j=i,min(i+5,maxlab)
               jcount=jcount+1
               TotEl(jcount)=PrEl(1,j)
            End Do
         Else
            Do j=1,6
               TotEl(j)=Zero
            End Do
            Write (outlab,'(1x,a,6a16)') 'Irrep  Orb   Occupation',
     &                                (labs(j),j=i,min(i+5,maxlab))
            Write(6,'(a)') outlab(:mylen(outlab))
            lbb4='    '
            Write (outlab,'(33a4)') (lbb4,j=1,33)
            lbb4='----'
            lab ='----------------'
            Write (outlab,'(6a4,6a16)') (lbb4,j=1,6),
     &                              (lab,j=i,min(i+5,maxlab))
            Write(6,'(1x,a)') outlab(:mylen(outlab))
            icount=0
            Do ii=0,nIrrep-1
               Do jj=1,nBas(ii)
                  icount=icount+1
                  jcount=0
                  Do j=i,min(i+5,maxlab)
                     jcount=jcount+1
                      TotEl(jcount)=TotEl(jcount)+PrEl(icount,j)
                  End Do
                  If (Occ(icount).gt.ThrSV)
     &               Write (6,Format1) ii+1,jj,Occ(icount),
     &                         (sig*PrEl(icount,j),j=i,min(i+5,maxlab))
               End Do
            End Do
            Write(6,'(1x,a)') outlab(:mylen(outlab))
         End if
*

         Do j=i,min(i+5,maxlab)
            PrTot(j)=sig*TotEl(j-i+1)+PrNu(j)
         End Do
         If (iPL.ge.3.or.(.Not.Short.and.iPL.eq.2)) Then
           Write(6,Format2)    'Total electronic       ',
     &              (sig*TotEl(j-i+1),j=i,min(i+5,maxlab))
           Write(6,Format2)    'Total nuclear          ',
     &              (PrNu(j),j=i,min(i+5,maxlab))
           Write(6,Format3)    'Total                  ',
     &              (PrTot(j),j=i,min(i+5,maxlab))
         Else
           lab5='     '
           If (iPoint.gt.0.and.i.eq.1) Write(lab5,'(I5)') iPoint
           Write(6,Format4) lab5,(PrTot(j),j=i,min(i+5,maxlab))
         End If
      End Do
*
      Return
      End
