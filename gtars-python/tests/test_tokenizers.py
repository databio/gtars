import os
from pathlib import Path

import pytest

from gtars.tokenizers import Tokenizer
from gtars.models import Region

TEST_DATA_DIR = os.path.abspath(
    os.path.join(Path(__file__).parents[3], "gtars/tests/data/")
)


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
    assert tokenizer.vocab_size == 32  # 25 + 7 special tokens


def test_tokenizer_creation_from_bed(tokenizer_peaks_bed: str):
    tokenizer = Tokenizer.from_bed(tokenizer_peaks_bed)
    assert tokenizer is not None
    assert tokenizer.vocab_size == 32  # 25 + 7 special tokens


def test_tokenizer_creation_from_bed_gz(tokenizer_peaks_bed_gz: str):
    tokenizer = Tokenizer.from_bed(tokenizer_peaks_bed_gz)
    assert tokenizer is not None
    assert tokenizer.vocab_size == 32  # 25 + 7 special tokens


@pytest.mark.parametrize(
    "path",
    [
        os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.bed"),
        os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.bed.gz"),
        os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml"),
        os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_ordered.toml"),
        os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_custom_specials.toml"),
    ],
)
def test_tokenizer_creation_auto_all(path: str):
    tokenizer = Tokenizer(path)
    assert tokenizer is not None
    assert tokenizer.vocab_size == 32  # 25 + 7 special tokens


def test_tokenizer_bad_tokenizer_type():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer_bad_ttype.toml")
    with pytest.raises(Exception):
        Tokenizer.from_config(cfg_path)


def test_tokenizer_custom_special_tokens():
    cfg_path = os.path.join(
        TEST_DATA_DIR, "tokenizers", "tokenizer_custom_specials.toml"
    )
    tokenizer = Tokenizer.from_config(cfg_path)
    assert tokenizer is not None
    assert tokenizer.vocab_size == 32  # 25 + 7 special tokens

    # check that unk was overridden
    unk_token = tokenizer.unk_token
    assert unk_token == "<UNKNOWN>"

    # check that pad didn't change
    pad_token = tokenizer.pad_token
    assert pad_token == "<pad>"


def test_tokenize_single_region_not_overlapping():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [Region("chr1", 50, 150, None)]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized) == 1
    assert tokenizer.convert_tokens_to_ids(tokenized[0]) == 25


def test_tokenize_unk_chrom():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [Region("chr999", 50, 150, None)]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized) == 1
    assert tokenizer.convert_tokens_to_ids(tokenized[0]) == 25


def test_tokenize_on_two_chroms():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [
        Region("chr1", 151399441, 151399547, None),
        Region("chr2", 203871220, 203871381, None),
    ]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized) == 2

    assert tokenized[0] == "chr1:151399431-151399527"
    assert tokenizer.convert_tokens_to_ids(tokenized[0]) == 6

    assert tokenized[1] == "chr2:203871200-203871375"
    assert tokenizer.convert_tokens_to_ids(tokenized[1]) == 7


def test_tokenize_with_multi_overlap():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "tokenizer.toml")
    tokenizer = Tokenizer.from_config(cfg_path)
    regions = [Region("chr2", 203871346, 203871616, None)]
    tokenized = tokenizer.tokenize(regions)
    assert tokenized is not None
    assert len(tokenized) == 2

    assert tokenized[0] == "chr2:203871200-203871375"
    assert tokenizer.convert_tokens_to_ids(tokenized[0]) == 7

    assert tokenized[1] == "chr2:203871387-203871588"
    assert tokenizer.convert_tokens_to_ids(tokenized[1]) == 8


def test_get_vocab():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = Tokenizer(cfg_path)

    vocab = tokenizer.get_vocab()
    assert vocab is not None
    assert len(vocab) == 32  # 25 + 7 special tokens
    assert all([isinstance(k, str) and isinstance(v, int) for k, v in vocab.items()])


def test_special_tokens_map():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = Tokenizer(cfg_path)

    special_tokens_map = tokenizer.special_tokens_map
    assert special_tokens_map is not None
    assert isinstance(special_tokens_map, dict)
    assert all(
        [
            isinstance(k, str) and isinstance(v, str)
            for k, v in special_tokens_map.items()
        ]
    )


def test_encode_tokens():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = Tokenizer(cfg_path)

    encoded = tokenizer.encode("chr9:3526071-3526165")
    assert encoded is not None
    assert isinstance(encoded, list)
    assert len(encoded) == 1
    assert encoded == [11]


def test_decode_tokens():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = Tokenizer(cfg_path)

    decoded = tokenizer.decode([11])
    assert decoded is not None
    assert isinstance(decoded, list)
    assert decoded == ["chr9:3526071-3526165"]


# @pytest.mark.skip(reason="Needs to be fixed")
def test_special_tokens_mask():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = Tokenizer(cfg_path)

    tokens = ["<pad>", "chr9:3526071-3526165", "<unk>"]
    special_tokens_mask = tokenizer.get_special_tokens_mask(tokens)

    assert special_tokens_mask is not None
    assert isinstance(special_tokens_mask, list)
    assert len(special_tokens_mask) == 3
    assert special_tokens_mask == [
        1,
        0,
        1,
    ]  # Assuming <pad> and <unk> are special tokens


def test_tokenizer_call_magic_method():
    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = Tokenizer(cfg_path)

    rs = [Region("chr9", 3526178, 3526249, None)]
    encoded = tokenizer(rs)

    assert encoded["input_ids"] == [10]
    assert encoded["attention_mask"] == [1]
    with pytest.raises(Exception):
        encoded["token_type_ids"]


def test_tokenizer_is_subclassable():
    class MoreTokenizer(Tokenizer):
        def __new__(cls, *args, **kwargs):
            return super().__new__(cls, *args, **kwargs)

        def __init__(self, *args, **kwargs):
            super().__init__()

        def add_two(self, x, y):
            return x + y

        @property
        def value(self):
            return self._value

        @value.setter
        def value(self, val):
            self._value = val

    cfg_path = os.path.join(TEST_DATA_DIR, "tokenizers", "peaks.scored.bed")
    tokenizer = MoreTokenizer(cfg_path)
    tokenizer.value = 5

    rs = [Region("chr9", 3526178, 3526249, None)]
    encoded = tokenizer(rs)

    assert encoded["input_ids"] == [10]
    assert encoded["attention_mask"] == [1]
    assert tokenizer.add_two(2, 3) == 5
    assert tokenizer.value == 5
