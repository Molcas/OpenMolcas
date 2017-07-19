#ifndef _HAVE_EXTRA_

      Subroutine ProcessXYZ(BasisSet, KeepBasis, KeepGroup, iBSSE,
     &                      SymThr, isHold, ScaleFactor, HyperParSet,
     &                      isXfield)
      Logical :: BasisSet, HyperParSet
      Character(Len=*) :: KeepBasis, KeepGroup
      Integer :: iBSSE, isHold, isXfield
      Real*8 :: SymThr, ScaleFactor
      Call Unused_Logical(BasisSet)
      Call Unused_Character(KeepBasis)
      Call Unused_Character(KeepGroup)
      Call Unused_Integer(iBSSE)
      Call Unused_Real(SymThr)
      Call Unused_Integer(isHold)
      Call Unused_Real(ScaleFactor)
      Call Unused_Logical(HyperParSet)
      Call Unused_Integer(isXfield)
      End Subroutine ProcessXYZ

#endif
