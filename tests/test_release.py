from __future__ import annotations

import configparser
import json
from pathlib import Path
from xml.etree import ElementTree

from StORF_Reporter import __version__
from StORF_Reporter.io import read_fasta, read_gff

ROOT = Path(__file__).parents[1]


def test_release_versions_and_remover_entrypoint_are_consistent() -> None:
    config = configparser.ConfigParser()
    config.read(ROOT / "setup.cfg")
    assert config["metadata"]["version"] == "attr: StORF_Reporter.constants.__version__"
    assert __version__ == "2.0.0"
    citation = (ROOT / "CITATION.cff").read_text()
    assert f'version: "{__version__}"' in citation
    assert json.loads((ROOT / "metadata/biotools.json").read_text())["version"] == [__version__]
    ro_crate = json.loads((ROOT / "metadata/ro-crate-metadata.json").read_text())
    software = next(item for item in ro_crate["@graph"] if item.get("@id") == "./")
    assert software["version"] == __version__
    entrypoints = config["options.entry_points"]["console_scripts"]
    assert "StORF-Remover = StORF_Reporter.StORF_Remover:main" in entrypoints
    assert "storf = StORF_Reporter.cli:main" in entrypoints


def test_supported_source_does_not_use_eval_or_hardcoded_home_paths() -> None:
    source = "\n".join(path.read_text() for path in (ROOT / "src").rglob("*.py"))
    assert "type=eval" not in source.replace(" ", "")
    assert "/home/" not in source


def test_quickstart_is_nonempty_and_coordinate_consistent() -> None:
    fasta = read_fasta(ROOT / "examples/quickstart/genome.fna")
    gff = read_gff(ROOT / "examples/quickstart/genome.gff3")
    assert len(fasta["ctg"]) == 229
    assert gff.features[0].interval.end <= len(fasta["ctg"])


def test_machine_readable_metadata_and_galaxy_xml_parse() -> None:
    json.loads((ROOT / "metadata/biotools.json").read_text())
    json.loads((ROOT / "metadata/ro-crate-metadata.json").read_text())
    ElementTree.parse(ROOT / "packaging/galaxy/storf_augment.xml")


def test_release_record_contains_scientific_warning_and_gates() -> None:
    release = (ROOT / "RELEASE_2.0.0.md").read_text()
    assert "not a" in release.lower() and "probability" in release.lower()
    assert "79 passed" in release
    assert "Publication benchmark" in release
    assert "Circular replicons" in release
