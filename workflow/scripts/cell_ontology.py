"""
MatchACell - minimal Cell Ontology (CL) reader for the consensus step.

Reads an OBO file (cl-basic.obo) and keeps, for every non-obsolete CL term,
its name, its EXACT synonyms and its is_a parents. Only is_a is followed:
"T cell" is an ancestor of "CD4-positive, alpha-beta T cell", but a pro-B cell
is not a descendant of "B cell" (CL links the two by development only).
"""

from __future__ import annotations

import re
from functools import lru_cache

_SYNONYM = re.compile(r'"(.*)" (EXACT|RELATED|BROAD|NARROW)')


def normalize(label) -> str:
    """Case-, spacing- and plural-insensitive key for matching labels."""
    s = str(label).lower().strip()
    s = re.sub(r"[\s_]+", " ", s)
    s = re.sub(r"\bcells\b", "cell", s)
    # singular last word ("monocytes", "T lymphocytes"), on labels and CL names
    # alike; keeps -ss/-us/-is endings ("plexus", "mass")
    return re.sub(r"(\w{3,}[^sui])s$", r"\1", s)


class CellOntology:
    def __init__(self, path):
        self.names, self.parents, synonyms = {}, {}, {}
        self._by_label = {}  # names only
        self.version = ""
        term = None
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("data-version: ") and not self.version:
                    self.version = line.split(": ", 1)[1]
                if line.startswith("["):
                    self._add(term, synonyms)
                    term = {"is_a": [], "syn": []} if line == "[Term]" else None
                    continue
                if term is None or ": " not in line:
                    continue
                key, value = line.split(": ", 1)
                if key == "id":
                    term["id"] = value
                elif key == "name":
                    term["name"] = value
                elif key == "is_a":
                    term["is_a"].append(value.split(" ")[0])
                elif key == "synonym":
                    m = _SYNONYM.match(value)
                    if m and m.group(2) == "EXACT":
                        term["syn"].append(m.group(1))
                elif key == "is_obsolete" and value == "true":
                    term["obsolete"] = True
        self._add(term, synonyms)
        self.parents = {t: [p for p in ps if p in self.names] for t, ps in self.parents.items()}
        self._by_name = synonyms
        self.ancestors = lru_cache(maxsize=None)(self._ancestors)

    def _add(self, term, synonyms):
        if not term or not term.get("id", "").startswith("CL:") or term.get("obsolete"):
            return
        tid = term["id"]
        self.names[tid] = term.get("name", tid)
        self.parents[tid] = term["is_a"]
        self._by_label.setdefault(normalize(self.names[tid]), set()).add(tid)
        for s in [term.get("name", "")] + term["syn"]:
            synonyms.setdefault(normalize(s), set()).add(tid)

    def __contains__(self, tid):
        return tid in self.names

    def _ancestors(self, tid) -> frozenset:
        """`tid` and every is_a ancestor of it."""
        out = {tid}
        for p in self.parents.get(tid, ()):
            out |= self.ancestors(p)
        return frozenset(out)

    def specificity(self, tid) -> int:
        """Number of ancestors: higher means a more specific term."""
        return len(self.ancestors(tid))

    def lookup(self, label):
        """CL id whose name or EXACT synonym matches `label`, else None
        (also None when the match is ambiguous). Abbreviations of up to four
        characters match names only: as synonyms they are ambiguous ("DC2" is
        an old name of plasmacytoid dendritic cells)."""
        key = normalize(label)
        index = self._by_label if len(key) <= 4 else self._by_name
        for k in (key, key + " cell"):
            ids = index.get(k)
            if ids and len(ids) == 1:
                return next(iter(ids))
        return None
