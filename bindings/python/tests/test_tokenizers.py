import os

import pytest

from gtars.tokenizers import Tokenizer

TEST_DATA_DIR = os.path.abspath('../../gtars/tests/data/')

@pytest.fixture
def tokenizer_config() -> str:
    return os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")

def test_tokenizer_initialization(tokenizer_config: str):
    tokenizer = Tokenizer.from_config(tokenizer_config)

    assert tokenizer is not None
    assert tokenizer.get_vocab_size() == 32 # 25 + 7 special tokens