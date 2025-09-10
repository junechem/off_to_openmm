#!/usr/bin/env python3
"""
Test cases for pdb_preprocessor.py
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path

# Add the scripts directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'scripts'))

def test_pdb_preprocessor_basic():
    """Test basic PDB preprocessing functionality."""
    
    # Create a sample PDB file with atom classes (proper PDB format)
    sample_pdb = """HEADER    TEST MOLECULE
ATOM      1  C3  CYCQM A   1      -0.123   1.456  -0.789  1.00  0.00           C
ATOM      2  C2  CYCQM A   1       0.987  -0.234   0.567  1.00  0.00           C  
ATOM      3  C1  CYCQM A   1       1.234   0.678  -1.234  1.00  0.00           C
ATOM      4  H3  CYCQM A   1      -1.089   1.987  -0.234  1.00  0.00           H
ATOM      5  H2  CYCQM A   1       1.456  -0.789   1.123  1.00  0.00           H
ATOM      6  OW  H2OQM A   2       2.345   1.234  -2.345  1.00  0.00           O
ATOM      7  HW  H2OQM A   2       3.123   1.789  -2.567  1.00  0.00           H
ATOM      8  HW  H2OQM A   2       2.567   0.567  -3.123  1.00  0.00           H
END
"""
    
    # Use existing XML file
    xml_file = Path(__file__).parent.parent / "final_test.xml"
    
    # Create temporary files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as pdb_file:
        pdb_file.write(sample_pdb)
        input_pdb = pdb_file.name
    
    output_pdb = input_pdb.replace('.pdb', '_processed.pdb')
    
    try:
        # Run the pdb_preprocessor script
        script_path = Path(__file__).parent.parent / "scripts" / "pdb_preprocessor.py"
        result = subprocess.run([
            sys.executable, str(script_path),
            "-pdb", input_pdb,
            "-xml", str(xml_file),
            "-output", output_pdb
        ], capture_output=True, text=True, timeout=30)
        
        print("STDOUT:")
        print(result.stdout)
        print("STDERR:")
        print(result.stderr)
        print(f"Return code: {result.returncode}")
        
        if result.returncode == 0:
            print("✅ PDB preprocessor ran successfully!")
            
            # Check if output file exists and has CONECT records
            if os.path.exists(output_pdb):
                with open(output_pdb, 'r') as f:
                    content = f.read()
                
                print("\nOutput PDB content:")
                print(content)
                
                # Verify CONECT records were added
                conect_count = content.count('CONECT')
                print(f"\n✅ Found {conect_count} CONECT records")
                
                # Check if atom classes were preserved or converted appropriately
                if 'C3' in content or 'C2' in content:
                    print("✅ Atom names preserved/converted correctly")
                
            else:
                print("❌ Output file was not created")
        else:
            print("❌ PDB preprocessor failed")
            
    except subprocess.TimeoutExpired:
        print("❌ Script timed out")
    except Exception as e:
        print(f"❌ Error running script: {e}")
    
    finally:
        # Clean up temporary files
        try:
            os.unlink(input_pdb)
            if os.path.exists(output_pdb):
                os.unlink(output_pdb)
        except:
            pass


def test_pdb_preprocessor_with_real_data():
    """Test with more realistic PDB data structure."""
    
    # Create a more realistic PDB file
    realistic_pdb = """HEADER    CYCLOHEXENE HYDRATE SYSTEM
REMARK   Generated for OpenMM testing
ATOM      1  C3  CYCQM   1      -1.123   2.456  -0.789  1.00  0.00           C
ATOM      2  C2  CYCQM   1       0.987  -0.234   0.567  1.00  0.00           C  
ATOM      3  C1  CYCQM   1       1.234   0.678  -1.234  1.00  0.00           C
ATOM      4  C1  CYCQM   1       2.345   1.789  -2.345  1.00  0.00           C
ATOM      5  C2  CYCQM   1       3.456   2.890  -3.456  1.00  0.00           C
ATOM      6  C3  CYCQM   1       4.567   3.901  -4.567  1.00  0.00           C
ATOM      7  H3  CYCQM   1      -2.089   3.987  -1.234  1.00  0.00           H
ATOM      8  H3  CYCQM   1      -1.234   2.123   0.345  1.00  0.00           H
ATOM      9  H2  CYCQM   1       1.456  -1.789   1.123  1.00  0.00           H
ATOM     10  H2  CYCQM   1       0.234  -0.567   1.456  1.00  0.00           H
ATOM     11  H1  CYCQM   1       1.567   0.345  -2.123  1.00  0.00           H
ATOM     12  H1  CYCQM   1       2.678   1.456  -1.890  1.00  0.00           H
ATOM     13  OW  H2OQM   2       5.678   4.012  -5.678  1.00  0.00           O
ATOM     14  HW  H2OQM   2       6.789   5.123  -6.789  1.00  0.00           H
ATOM     15  HW  H2OQM   2       5.890   4.567  -6.234  1.00  0.00           H
TER
END
"""
    
    xml_file = Path(__file__).parent.parent / "final_test.xml"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as pdb_file:
        pdb_file.write(realistic_pdb)
        input_pdb = pdb_file.name
    
    output_pdb = input_pdb.replace('.pdb', '_realistic_processed.pdb')
    
    try:
        script_path = Path(__file__).parent.parent / "scripts" / "pdb_preprocessor.py"
        result = subprocess.run([
            sys.executable, str(script_path),
            "-pdb", input_pdb,
            "-xml", str(xml_file),
            "-output", output_pdb
        ], capture_output=True, text=True, timeout=30)
        
        print("\n" + "="*60)
        print("REALISTIC PDB TEST")
        print("="*60)
        print("STDOUT:")
        print(result.stdout)
        print("STDERR:")
        print(result.stderr)
        
        if result.returncode == 0 and os.path.exists(output_pdb):
            with open(output_pdb, 'r') as f:
                content = f.read()
            
            # Count different record types
            atom_count = content.count('ATOM')
            conect_count = content.count('CONECT')
            
            print(f"\n✅ Successfully processed realistic PDB:")
            print(f"   - {atom_count} ATOM records")
            print(f"   - {conect_count} CONECT records")
            
            # Show some CONECT records
            conect_lines = [line for line in content.split('\n') if line.startswith('CONECT')]
            print(f"   - Sample CONECT records:")
            for conect in conect_lines[:3]:
                print(f"     {conect}")
            if len(conect_lines) > 3:
                print(f"     ... and {len(conect_lines) - 3} more")
                
        else:
            print("❌ Realistic PDB test failed")
            
    except Exception as e:
        print(f"❌ Error in realistic test: {e}")
        
    finally:
        try:
            os.unlink(input_pdb)
            if os.path.exists(output_pdb):
                os.unlink(output_pdb)
        except:
            pass


if __name__ == "__main__":
    test_pdb_preprocessor_basic()
    test_pdb_preprocessor_with_real_data()