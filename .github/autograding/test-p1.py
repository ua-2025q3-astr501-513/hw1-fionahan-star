# This is a pytest script to test the package import for Part 1.
#
# It ensures that the module structure and dependencies
# (pyproject.toml, etc.) are set up correctly so that the required
# submodules can be imported.

from pathlib import Path

def test_import():
    """Test whether the module structure is correctly configured for
    imports.

    This test checks if the package `hw1` is properly installed and if
    its submodules (p2, p3, p4) can be successfully imported.

    If `pyproject.toml` and other configuration files are correctly
    set up, the following import should work without errors.

    """
    try:
        from hw1 import p2, p3, p4  # import required submodules
    except ImportError as e:
        assert False, f"Module import failed: {e}"

    if not Path('LICENSE').is_file():
        assert False, "'LICENSE' does not exist"
