import pytest
import os
import subprocess
import tempfile
import shutil
from pathlib import Path


class TestOffToOpenMMIntegration:
    """Integration test suite for off_to_openmm.py with real .off files."""
    
    @pytest.fixture
    def test_data_dir(self):
        """Fixture to provide path to test data directory."""
        return Path(__file__).parent / "test_data"
    
    @pytest.fixture
    def temp_output_dir(self):
        """Fixture to provide temporary output directory."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    def test_intra_off_conversion(self, test_data_dir, temp_output_dir):
        """Test conversion of intra.off file to OpenMM XML."""
        # Copy intra.off to test data directory for testing
        input_file = test_data_dir / "intra.off"
        output_file = Path(temp_output_dir) / "test_output.xml"
        
        # Skip if test file doesn't exist yet
        if not input_file.exists():
            pytest.skip("intra.off test file not available yet")
        
        # Run the conversion script
        result = subprocess.run([
            "python", "off_to_openmm.py",
            str(input_file),
            str(output_file)
        ], capture_output=True, text=True)
        
        # Check that script ran successfully
        assert result.returncode == 0, f"Script failed with error: {result.stderr}"
        
        # Check that output file was created
        assert output_file.exists(), "Output XML file was not created"
        
        # Basic XML structure validation
        with open(output_file, 'r') as f:
            content = f.read()
            assert "<AtomTypes>" in content, "AtomTypes section missing"
            assert "<Residues>" in content, "Residues section missing"
            assert "<HarmonicBondForce>" in content, "HarmonicBondForce section missing"
    
    def test_multiple_off_files(self, test_data_dir, temp_output_dir):
        """Test conversion with multiple different .off files."""
        test_files = list(test_data_dir.glob("*.off"))
        
        if len(test_files) == 0:
            pytest.skip("No .off test files available yet")
        
        for off_file in test_files:
            output_file = Path(temp_output_dir) / f"{off_file.stem}_output.xml"
            
            # Run conversion
            result = subprocess.run([
                "python", "off_to_openmm.py", 
                str(off_file),
                str(output_file)
            ], capture_output=True, text=True)
            
            assert result.returncode == 0, f"Conversion failed for {off_file.name}: {result.stderr}"
            assert output_file.exists(), f"Output not created for {off_file.name}"
    
    def test_command_line_flags(self, test_data_dir, temp_output_dir):
        """Test command line flags like -molnames and -charges."""
        input_file = test_data_dir / "intra.off"
        output_file = Path(temp_output_dir) / "test_flags.xml"
        
        if not input_file.exists():
            pytest.skip("intra.off test file not available yet")
        
        # Test -molnames flag
        result = subprocess.run([
            "python", "off_to_openmm.py",
            str(input_file),
            str(output_file),
            "-molnames", "UNK"
        ], capture_output=True, text=True)
        
        assert result.returncode == 0, f"Script failed with -molnames flag: {result.stderr}"


class TestOffToOpenMMUnit:
    """Unit tests for individual components of off_to_openmm.py."""
    
    def test_off_file_parsing(self):
        """Test parsing of .off file format."""
        # TODO: Implement once off_to_openmm.py is created
        pass
    
    def test_xml_generation(self):
        """Test OpenMM XML format generation."""
        # TODO: Implement once off_to_openmm.py is created
        pass