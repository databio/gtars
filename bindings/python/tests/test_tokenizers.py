import os
from pathlib import Path

import pytest

from gtars.tokenizers import Tokenizer
from gtars.models import Region

TEST_DATA_DIR = os.path.abspath(os.path.join(Path(__file__).parents[3], 'gtars/tests/data/'))

@pytest.fixture
def tokenizer_config() -> str:
    return os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")

@pytest.fixture
def tokenizer_bad_config() -> str:
    return os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_bad.toml")

@pytest.fixture
def tokenizer_custom_specials_config() -> str:
    return os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_custom_specials.toml")

@pytest.fixture
def tokenizer_bad_ttype_config() -> str:
    return os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_bad_ttype.toml")

@pytest.fixture
def tokenizer_peaks_bed() -> str:
    return os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.bed")

@pytest.fixture
def tokenizer_peaks_bed_gz() -> str:
    return os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.bed.gz")

def test_tokenizer_initialization(tokenizer_config: str):
    print(tokenizer_config)
    tokenizer = Tokenizer.from_config(tokenizer_config)

    assert tokenizer is not None
    assert tokenizer.get_vocab_size() == 32 # 25 + 7 special tokens
    
def test_tokenizer_creation_from_bed(tokenizer_peaks_bed: str):
    tokenizer = Tokenizer.from_bed(tokenizer_peaks_bed)
    assert tokenizer is not None
    assert tokenizer.get_vocab_size() == 32 # 25 + 7 special tokens

def test_tokenizer_creation_from_bed_gz(tokenizer_peaks_bed_gz: str):
    tokenizer = Tokenizer.from_bed(tokenizer_peaks_bed_gz)
    assert tokenizer is not None
    assert tokenizer.get_vocab_size() == 32 # 25 + 7 special tokens

@pytest.mark.parametrize("path", [
    os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.bed"),
    os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.bed.gz"),
    os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml"),
    os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_ordered.toml"),
    os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_custom_specials.toml")
])  
def test_tokenizer_creation_auto_all(path: str):
    tokenizer = Tokenizer(path)
    assert tokenizer is not None
    assert tokenizer.get_vocab_size() == 32 # 25 + 7 special tokens  
    
def test_tokenizer_bad_tokenizer_type():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_bad_ttype.toml")
    with pytest.raises(Exception):
        Tokenizer.from_config(cfg_path)

def test_tokenizer_custom_special_tokens():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_custom_specials.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    assert tokenizer is not None
    assert tokenizer.get_vocab_size() == 32 # 25 + 7 special tokens

    # check that unk was overridden
    unk_token = tokenizer.get_unk_token()
    assert unk_token.chr == "chrUNKNOWN"
    assert unk_token.start == 100
    assert unk_token.end == 200

    # check that pad didn't change
    pad_token = tokenizer.get_pad_token()
    assert pad_token.chr == "chrPAD"
    assert pad_token.start == 0
    assert pad_token.end == 0

def test_tokenize_single_region_not_overlapping():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [Region("chr1", 50, 150)]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized.ids) == 1
    assert tokenized.ids[0] == tokenizer.get_unk_token_id()
    assert tokenizer.convert_id_to_token(tokenized.ids[0]).chr == "chrUNK"

def test_tokenize_unk_chrom():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [Region("chr999", 50, 150)]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized.ids) == 1

def test_tokenize_on_two_chroms():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [
        Region("chr1", 151399441, 151399547),
        Region("chr2", 203871220, 203871381)
    ]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized) == 2

    tokenized_regions = tokenized.to_regions()
    assert tokenized_regions[0].chr == "chr1"
    assert tokenized_regions[0].start == 151399431
    assert tokenized_regions[0].end == 151399527
    assert tokenizer.convert_token_to_id(tokenized_regions[0]) == 6

    assert tokenized_regions[1].chr == "chr2"
    assert tokenized_regions[1].start == 203871200
    assert tokenized_regions[1].end == 203871375
    assert tokenizer.convert_token_to_id(tokenized_regions[1]) == 7

def test_tokenize_with_multi_overlap():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [Region("chr2", 203871346, 203871616)]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized) == 2

    tokenized_regions = tokenized.to_regions()
    assert tokenized_regions[0].chr == "chr2"
    assert tokenized_regions[0].start == 203871200
    assert tokenized_regions[0].end == 203871375
    assert tokenizer.convert_token_to_id(tokenized_regions[0]) == 7

    assert tokenized_regions[1].chr == "chr2"
    assert tokenized_regions[1].start == 203871387
    assert tokenized_regions[1].end == 203871588
    assert tokenizer.convert_token_to_id(tokenized_regions[1]) == 8

def test_tokenize_with_order():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = Tokenizer(cfg_path)
    regions = [
        Region("chr9", 3526178, 3526249),
        Region("chr9", 3526051, 3526145)
    ]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized) == 2

    tokenized_regions = tokenized.to_regions()
    assert tokenized_regions[0].chr == "chr9"
    assert tokenized_regions[0].start == 3526071
    assert tokenized_regions[0].end == 3526165
    assert tokenizer.convert_token_to_id(tokenized_regions[0]) == 11

    assert tokenized_regions[1].chr == "chr9"
    assert tokenized_regions[1].start == 3526183
    assert tokenized_regions[1].end == 3526269
    assert tokenizer.convert_token_to_id(tokenized_regions[1]) == 18
