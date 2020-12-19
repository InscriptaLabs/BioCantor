import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir() -> Path:
    return Path(__file__).parent / "data"
