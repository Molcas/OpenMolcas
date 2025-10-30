#!/usr/bin/env python3
"""
Quick test to verify GA4GUGA package installation
Run this after: pip install -e .
"""

import sys

def test_imports():
    """Test that all modules can be imported"""
    print("Testing GA4GUGA package installation...\n")
    
    tests_passed = 0
    tests_failed = 0
    
    # Test GA_mod
    try:
        import GA_mod
        print("✓ GA_mod imported successfully")
        print(f"  Location: {GA_mod.__file__}")
        tests_passed += 1
    except ImportError as e:
        print(f"✗ Failed to import GA_mod: {e}")
        tests_failed += 1
    
    # Test GA_mod submodules
    try:
        from GA_mod import run_GA
        from GA_mod import crossover
        from GA_mod import population
        from GA_mod.measure_fitness import FitnessFunction
        print("✓ GA_mod submodules imported successfully")
        tests_passed += 1
    except ImportError as e:
        print(f"✗ Failed to import GA_mod submodules: {e}")
        tests_failed += 1
    
    # Test FCIDUMP_tools
    try:
        import FCIDUMP_tools
        print("✓ FCIDUMP_tools imported successfully")
        print(f"  Location: {FCIDUMP_tools.__file__}")
        tests_passed += 1
    except ImportError as e:
        print(f"✗ Failed to import FCIDUMP_tools: {e}")
        tests_failed += 1
    
    # Test FCIDUMP_tools submodules
    try:
        from FCIDUMP_tools import IntegralClass
        from FCIDUMP_tools import xyz2heisenberg
        print("✓ FCIDUMP_tools submodules imported successfully")
        tests_passed += 1
    except ImportError as e:
        print(f"✗ Failed to import FCIDUMP_tools submodules: {e}")
        tests_failed += 1
    
    # Test dependencies
    try:
        import numpy
        import pandas
        print("✓ Required dependencies (numpy, pandas) available")
        tests_passed += 1
    except ImportError as e:
        print(f"✗ Missing required dependencies: {e}")
        tests_failed += 1
    
    print(f"\n{'='*50}")
    print(f"Results: {tests_passed} passed, {tests_failed} failed")
    print(f"{'='*50}")
    
    if tests_failed == 0:
        print("\n✓ All tests passed! Package is correctly installed.")
        print("\nYou can now:")
        print("  1. Remove PYTHONPATH from ~/.bash_profile")
        print("  2. Use 'from GA_mod import ...' in your scripts")
        return 0
    else:
        print("\n✗ Some tests failed. Try:")
        print("  cd /Users/song/CODE/GA4GUGA")
        print("  pip install -e .")
        return 1

if __name__ == "__main__":
    sys.exit(test_imports())
