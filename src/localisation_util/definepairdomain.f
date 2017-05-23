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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine DefinePairDomain(irc,iPairDomain,iClass,Rmin,iDomain,
     &                            RThr,Coord,nAtom,nOcc,nRThr)
************************************************************
*
*   <DOC>
*     <Name>DefinePairDomain</Name>
*     <Syntax>Call DefinePairDomain(irc,iPairDomain,iClass,Rmin,iDomain,
*                                  RThr,Coord,nAtom,nOcc,nRThr)</Syntax>
*     <Arguments>
*       \Argument{irc}{Return code}{Integer}{out}
*       \Argument{iPairDomain}{Pair domain definition}
*                {Integer(0:nAtom,nOcc*(nOcc+1)/2)}{out}
*       \Argument{iClass}{Classification}{Integer(nOcc*(nOcc+1)/2)}{out}
*       \Argument{Rmin}{Minimum interatomic distance in each pair}
*                      {Real*8(nOcc*(nOcc+1)/2)}{out}
*       \Argument{iDomain}{Orbital domain definition}
*                         {Integer(0:nAtom,nOcc)}{in}
*       \Argument{RThr}{Thresholds for classification}{Real*8(nRThr)}
*                      {in}
*       \Argument{Coord}{Nuclear coordinates}{Real*8(3,nAtom)}{in}
*       \Argument{nAtom}{Number of atoms}{Integer}{in}
*       \Argument{nOcc}
*                {Number of orbitals for which domains are defined}
*                {Integer}{in}
*       \Argument{nRThr}{Number of thresholds for classification}
*                       {Integer}{in}
*     </Arguments>
*     <Purpose>Define pair domains</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Pair domains are defined by the union of individual orbital
*        domains, iDomain(0:nAtom,nOcc), see SubRoutine DefineDomain.
*        On exit, the contents of iPairDomain array are:
*        iPairDomain(0,ij): number of atoms in pair domain ij.
*        iPairDomain(n,ij): id of atom n (=1,2,...,iPairDomain(0,ij))
*                           in pair domain ij.
*        The contents of Rmin array are:
*        Rmin(ij): minimum distance between atoms in pair ij.
*        The contents of iClass array are (for i=1,2,3,...,nRThr):
*        iClass(ij)=i-1 implies Rmin(ij) <= RThr(i) and
*        iClass(ij)=nRThr implies Rmin(ij) > RThr(nRThr)
*        Note that if classification is not wanted, simply specify
*        nRThr=0 (in which case arrays iClass and RThr are not referenced
*        may be dummy arguments in the call to this routine).
*        The iDomain array should be set by SubRoutine DefineDomain
*        before calling this routine.
*        Return codes:
*        irc=0: all OK.
*        (at the moment, this is the only possible return code, but that
*        might change.)
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer iPairDomain(0:nAtom,*), iClass(*)
      Real*8  Rmin(*)
      Integer iDomain(0:nAtom,nOcc)
      Real*8  RThr(*), Coord(3,nAtom)
#include "WrkSpc.fh"

C     Set return code.
C     ----------------

      irc = 0
      If (nOcc .lt. 2) Return

C     Set pair domains as union of individual domains: [ij]=[i]U[j] for
C     i>=j.
C     -----------------------------------------------------------------

      nnOcc = nOcc*(nOcc+1)/2
      lT = (nAtom+1)*nnOcc
      Call iCopy(lT,0,0,iPairDomain,1)

      l_Union = nAtom*nOcc
      Call GetMem('Union','Allo','Inte',ip_Union,l_Union)

      Call iCopy(l_Union,0,0,iWork(ip_Union),1)
      kOff = ip_Union - 1
      Do i = 1,nOcc
         iOff = kOff + nAtom*(i-1)
         Do iA = 1,iDomain(0,i)
            iAtom = iDomain(iA,i)
            iWork(iOff+iAtom) = 1
         End Do
      End Do

      kOff = ip_Union - 1
      ij = 0
      Do j = 1,nOcc
         ij = ij + 1
         lD = iDomain(0,j) + 1
         Call iCopy(lD,iDomain(0,j),1,iPairDomain(0,ij),1) ! case i=j
         jOff = kOff + nAtom*(j-1)
         Do i = j+1,nOcc
            iOff = kOff + nAtom*(i-1)
            iCount = 0
            ij = ij + 1
            Do kAtom = 1,nAtom
               isThere = iWork(jOff+kAtom) + iWork(iOff+kAtom)
               If (isThere .gt. 0) Then
                  iCount = iCount + 1
                  iPairDomain(iCount,ij) = kAtom
               End If
            End Do
            iPairDomain(0,ij) = iCount
         End Do
      End Do

      Call GetMem('Union','Free','Inte',ip_Union,l_Union)

C     Set min. distance between any two atoms in pairs of domains.
C     ------------------------------------------------------------

      ij = 0
      Do j = 1,nOcc
         ij = ij + 1
         Rmin(ij) = 0.0d0 ! case i=j
         Do i = j+1,nOcc
            ij = ij + 1
            Rmin(ij) = 1.0d15
            Do jA = 1,iDomain(0,j)
               jAtom = iDomain(jA,j)
               Do iA = 1,iDomain(0,i)
                  iAtom = iDomain(iA,i)
                  R = sqrt((Coord(1,iAtom)-Coord(1,jAtom))**2
     &                    +(Coord(2,iAtom)-Coord(2,jAtom))**2
     &                    +(Coord(3,iAtom)-Coord(3,jAtom))**2)
                  Rmin(ij) = min(Rmin(ij),R)
               End Do
            End Do
         End Do
      End Do

C     Set classification for each pair domain.
C     Note that the caller must have ordered the thresholds according to
C     RThr(1) < RThr(2) < RThr(3) < ... < RThr(nRThr)
C     ------------------------------------------------------------------

      If (nRThr .gt. 0) Then
         Call iCopy(nnOcc,nRThr,0,iClass,1)
         Do ij = 1,nnOcc
            i = 0
            Do While (i .lt. nRThr)
               i = i + 1
               If (Rmin(ij) .le. RThr(i)) Then
                  iClass(ij) = i - 1
                  i = nRThr ! break while loop
               End If
            End Do
         End Do
      End If

      End
