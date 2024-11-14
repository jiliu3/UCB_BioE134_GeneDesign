import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

@pytest.fixture
def rbs_checker():
    r = InternalRBSChecker()
    r.initiate()
    return r

def test_rbs_sequences(rbs_checker):
    rbs_sequences = [
        ('AGGAGGTTTATG', True),
        ('ATTAGAAAACGGAGGAGGAGGATGTAGCAGCATTAT', False),
        ('ATTAGAAAACGGAGTAGGAGGATGTAGCAGCATTAT', False),
        ('ATTAGAAAACGGATTAGGAGGATGTAGCAGCATTAT', True),
        ('TGCCGTGGAGGTAAAATG', False),
        ('ATCTAGGAGTGCAAAATG', False),
        ('AAAAGGAGGAAAATG', True),
        ('AAAGGAGGAAAAATG', False),
        ('AAAAGGAGGAAAAAAATG', False),
        ('AAAGGAGGAAAAAAAATG', False)

    ]

    print("\n>> Testing rbs with expected results")
    for seq, expected in rbs_sequences:
        result, rbs = rbs_checker.run(seq)
        print(f"Sequence: {seq}, Expected: {expected}, Got: {result}, RBS: {rbs}")
        assert result == expected, f"Test failed for sequence: {seq}. Expected {expected} but got {result}."