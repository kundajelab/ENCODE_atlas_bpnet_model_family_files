"""encode_portal.py -- resolve + describe the downloadable ENCODE-portal objects for the
Kundaje-lab ENCODE Atlas BPNet-family models.

The index is the four family tables published at
  https://github.com/kundajelab/ENCODE_atlas_bpnet_model_family_files
(one TSV per family; row = a trained model by ENCODE accession, column = an output object type,
cell = the ENCODE file accession ENCFF... to download). This module fetches those tables from the
GitHub raw link (cached), falling back to the copy bundled in this package (data/encode_portal/)
when offline.

It underpins the acquisition ladder for the co-scientist -- cheapest rung first:
  1. derived objects (motifs / contributions / hits / signal) -> answer many questions with NO model,
  2. download + stage the trained model -> novel in-silico perturbations,
  3. train from scratch -> only when nothing is released.

Pure stdlib (urllib + csv); py3.8-safe.
"""
import csv
import json
import os
import time
import urllib.request

RAW_BASE = os.environ.get(
    "AGCB_ENCODE_ATLAS_RAW",
    "https://raw.githubusercontent.com/kundajelab/ENCODE_atlas_bpnet_model_family_files/main")
REPO_URL = "https://github.com/kundajelab/ENCODE_atlas_bpnet_model_family_files"

_FAMILIES = ["bpnet", "chrombpnet", "procapnet", "reporternet"]
_BUNDLED = os.path.join(os.path.dirname(__file__), "data", "encode_portal")
_CACHE = os.environ.get("AGCB_ENCODE_ATLAS_CACHE",
                        os.path.join(os.path.dirname(__file__), "data", "encode_portal", "_cache"))
_CACHE_TTL = int(os.environ.get("AGCB_ENCODE_ATLAS_TTL", str(24 * 3600)))

# canonical object name -> candidate table columns, in preference order (bigWig before tar for the
# streamable/region-queryable objects; the accessors pick the format they want).
OBJECTS = {
    "model":                ["models (tar)"],
    "motifs":               ["sequence motifs (tar)"],
    "motifs_report":        ["sequence motifs report (tar)"],
    "motif_hits":           ["sequence motifs instances (bigBed)", "sequence motifs instances (tar)"],
    "contrib_counts":       ["counts sequence contribution scores (bigWig)",
                             "counts sequence contribution scores (tar)"],
    "contrib_profile":      ["profile sequence contribution scores (bigWig)",
                             "profile sequence contribution scores (tar)"],
    "signal_predicted":     ["predicted signal profile (bigWig)", "normalized predicted signal profile (bigWig)",
                             "normalized predicted signal profile (plus strand) (bigWig)"],
    "signal_observed":      ["observed signal profile (bigWig)", "normalized observed signal profile (bigWig)",
                             "normalized observed signal profile (plus strand) (bigWig)"],
    "signal_biascorrected": ["normalized bias-corrected predicted signal profile (bigWig)",
                             "bias-corrected predicted signal profile (bigWig)"],
    "regions":              ["selected regions for predicted signal and sequence contribution scores (bed)"],
    "regions_train_test":   ["training and test regions (tar)"],
    "metrics":              ["model performance metrics (tar)"],
}
_META_COLS = 7   # experiment,type,assay,target,tissue,model_annotation,qc

_INDEX = None     # {accession: {"_family":.., "_meta":{...}, <column>: encff}}


def _read_table_text(fam):
    """Fetch one family table: GitHub raw (cached) -> bundled fallback. Returns TSV text."""
    url = RAW_BASE + "/" + fam + "_model_files.tsv"
    cache = os.path.join(_CACHE, fam + "_model_files.tsv")
    if os.path.exists(cache) and (time.time() - os.path.getmtime(cache)) < _CACHE_TTL:
        try:
            return open(cache, encoding="utf-8").read()
        except Exception:
            pass
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            text = r.read().decode("utf-8")
        try:
            os.makedirs(_CACHE, exist_ok=True)
            with open(cache, "w", encoding="utf-8") as fh:
                fh.write(text)
        except Exception:
            pass
        return text
    except Exception:
        # offline fallbacks: the bundled copy, then a sibling table (the public-repo layout, where
        # encode_portal.py sits next to the four *_model_files.tsv)
        for cand in (os.path.join(_BUNDLED, fam + "_model_files.tsv"),
                     os.path.join(os.path.dirname(os.path.abspath(__file__)), fam + "_model_files.tsv")):
            if os.path.exists(cand):
                return open(cand, encoding="utf-8").read()
        raise


def _load(force=False):
    global _INDEX
    if _INDEX is not None and not force:
        return _INDEX
    idx = {}
    for fam in _FAMILIES:
        try:
            rows = list(csv.reader(_read_table_text(fam).splitlines(), delimiter="\t"))
        except Exception:
            continue
        if not rows:
            continue
        hdr = rows[0]
        for r in rows[1:]:
            if not r or not r[0].strip():
                continue
            acc = r[0].strip()
            rec = idx.setdefault(acc, {"_family": fam, "_meta": {}})
            for i, col in enumerate(hdr):
                val = r[i].strip() if i < len(r) else ""
                if i < _META_COLS:
                    rec["_meta"][col] = val
                elif val:
                    rec[col] = val
    _INDEX = idx
    return idx


# ----------------------------------------------------------------- public API ---

def model_files(accession):
    """Everything released for one model: family, metadata, and the available object types with
    their ENCFF accessions. Returns None if the accession is not in the Atlas."""
    idx = _load()
    rec = idx.get(accession)
    if not rec:
        return None
    objects = {}
    for canon, cols in OBJECTS.items():
        for c in cols:
            if rec.get(c):
                objects[canon] = {"encff": rec[c], "column": c, "format": _fmt_of(c),
                                  "download_url": download_url(rec[c])}
                break
    return {"accession": accession, "family": rec["_family"], "meta": rec["_meta"],
            "objects": objects, "source": REPO_URL}


def resolve(accession, object_type, fmt=None):
    """The ENCFF accession for one (model, object_type). object_type is a canonical name from
    OBJECTS (e.g. 'motifs', 'contrib_counts', 'model'). Returns None if not released. If `fmt`
    is given ('bigWig'|'tar'|'bigBed'|'bed'), prefer a candidate column of that format."""
    idx = _load()
    rec = idx.get(accession)
    if not rec:
        return None
    cols = OBJECTS.get(object_type)
    if not cols:
        return None
    if fmt:
        for c in cols:
            if _fmt_of(c) == fmt and rec.get(c):
                return rec[c]
    for c in cols:
        if rec.get(c):
            return rec[c]
    return None


def download_url(encff):
    """The ENCODE portal download URL for a file accession (portal redirects to the real file)."""
    return "https://www.encodeproject.org/files/%s/@@download/" % encff


def file_metadata(encff, timeout=20):
    """ENCODE file metadata (format, size, md5, href) via the portal REST API."""
    url = "https://www.encodeproject.org/files/%s/?format=json" % encff
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        d = json.load(r)
    href = d.get("href")
    return {"accession": encff, "file_format": d.get("file_format"),
            "output_type": d.get("output_type"), "file_size": d.get("file_size"),
            "md5sum": d.get("md5sum"),
            "href": ("https://www.encodeproject.org" + href) if href else download_url(encff)}


def object_types():
    """The canonical object types this index can resolve, with a one-line description each."""
    return {
        "model": "trained model per fold (.h5 + SavedModel) -- run in-silico predictions/perturbations",
        "motifs": "TF-MoDISco motifs (CWM/PFM, MEME + h5) -- the learned motif vocabulary",
        "motifs_report": "TF-MoDISco HTML report",
        "motif_hits": "FiNeMo motif-instance calls (bigBed/tar) -- genome-wide occurrences",
        "contrib_counts": "per-base counts-head DeepSHAP attributions (bigWig/tar)",
        "contrib_profile": "per-base profile-head DeepSHAP attributions (bigWig/tar)",
        "signal_predicted": "predicted signal profile (bigWig)",
        "signal_observed": "measured signal profile (bigWig)",
        "signal_biascorrected": "ChromBPNet bias-corrected predicted signal (bigWig)",
        "regions": "regions the prediction/contribution files cover (bed)",
        "regions_train_test": "peaks/non-peaks per fold, train/valid/test split (tar)",
        "metrics": "per-fold Pearson/Spearman/AUROC/AUPRC/JSD (tar)",
    }


def _fmt_of(column):
    for f in ("bigWig", "bigBed", "tar", "bed"):
        if "(" + f + ")" in column or ("(" + f in column):
            return f
    return ""


# =============================================================================================
# General client -- usable standalone by ANY agent (pure stdlib). Download our models + their
# downstream products from the index, AND reach the broader ENCODE portal (which has far more
# than the BPNet-family models) via the public REST API.
# =============================================================================================

def download(accession, object_type, dest_dir=".", fmt=None, timeout=3600):
    """Download one released object (model / motifs / contributions / signal / hits / ...) to
    dest_dir. Returns the local file path. Raises if the object was not released for this model."""
    import urllib.request
    encff = resolve(accession, object_type, fmt=fmt)
    if not encff:
        raise ValueError("no %r object released for %s" % (object_type, accession))
    meta = {}
    try:
        meta = file_metadata(encff)
    except Exception:
        pass
    url = meta.get("href") or download_url(encff)
    ext = (meta.get("file_format") or "").lower()
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, ("%s.%s" % (encff, ext)) if ext else encff)
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=timeout) as r, open(tmp, "wb") as fh:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            fh.write(b)
    os.replace(tmp, dest)
    return dest


def search(family=None, assay=None, biosample=None, target=None, qc="passed", limit=200):
    """Filter OUR model index (the four family tables). Returns [{accession, family, assay, target,
    tissue, qc}]. For anything beyond our released models, use portal_search / portal_files."""
    idx = _load()
    out = []
    for acc, rec in idx.items():
        m, fam = rec["_meta"], rec["_family"]
        if family and fam != family.lower():
            continue
        if assay and assay.lower() not in (m.get("assay", "") or "").lower():
            continue
        if biosample and biosample.lower() not in (m.get("tissue", "") or "").lower():
            continue
        if target and target.lower() not in (m.get("target", "") or "").lower():
            continue
        if qc and (m.get("qc") or "") and qc.lower() not in (m.get("qc", "") or "").lower():
            continue
        out.append({"accession": acc, "family": fam,
                    "assay": m.get("assay"), "target": m.get("target"),
                    "tissue": m.get("tissue"), "qc": m.get("qc")})
        if len(out) >= limit:
            break
    return out


def portal_files(experiment, output_type=None, timeout=30):
    """BROADER than our models: every file released for ANY ENCODE experiment, via the public portal
    REST API. ENCODE hosts raw + processed data (reads, alignments, peaks, signal, etc.) for hundreds
    of thousands of experiments -- our tables index only the BPNet-family MODELS on top of a slice of
    them. Returns {experiment, assay, biosample, n_files, files:[{accession, output_type, file_format,
    assembly, download_url}]}."""
    import json as _json
    import urllib.request
    url = "https://www.encodeproject.org/experiments/%s/?format=json" % experiment
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        d = _json.load(r)
    files = []
    for f in d.get("files", []) or []:
        ot = f.get("output_type", "") or ""
        if output_type and output_type.lower() not in ot.lower():
            continue
        files.append({"accession": f.get("accession"), "output_type": ot,
                      "file_format": f.get("file_format"), "assembly": f.get("assembly"),
                      "download_url": download_url(f.get("accession", ""))})
    return {"experiment": experiment, "assay": d.get("assay_title"),
            "biosample": d.get("biosample_summary"), "n_files": len(files), "files": files}


def portal_search(assay_title=None, biosample=None, target=None, limit=25, timeout=30):
    """BROADER ENCODE experiment search (any experiment, not just our modelled ones) via the public
    REST API. Returns [{accession, assay, biosample, target}]."""
    import json as _json
    import urllib.parse
    import urllib.request
    q = [("type", "Experiment"), ("format", "json"), ("limit", str(limit))]
    if assay_title:
        q.append(("assay_title", assay_title))
    if biosample:
        q.append(("biosample_ontology.term_name", biosample))
    if target:
        q.append(("target.label", target))
    url = "https://www.encodeproject.org/search/?" + urllib.parse.urlencode(q)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        d = _json.load(r)
    return [{"accession": h.get("accession"), "assay": h.get("assay_title"),
             "biosample": h.get("biosample_summary"),
             "target": (h.get("target") or {}).get("label")} for h in d.get("@graph", []) or []]


# ------------------------------------------------- discovery (no ENCID needed) ---
# An agent that only knows "I want an accessibility model for liver" -- not the ENCODE accession --
# finds one via these. facets() shows what exists; find_models() free-text-ranks our index;
# find_models_by_organ() is ontology-aware (ENCODE's organ_slims maps 'liver' -> HepG2 / hepatocyte /
# liver tissue, which a plain tissue-string match would miss).

def facets():
    """What the index actually contains, so an agent can browse instead of guessing an accession:
    counts by family / assay / tissue / target. tissue+target sorted by frequency."""
    from collections import Counter
    idx = _load()
    fam, assay, tissue, target = Counter(), Counter(), Counter(), Counter()
    for rec in idx.values():
        m = rec["_meta"]
        fam[rec["_family"]] += 1
        for c, k in ((assay, "assay"), (tissue, "tissue"), (target, "target")):
            if m.get(k):
                c[m[k]] += 1
    return {"family": dict(fam), "assay": dict(assay.most_common()),
            "tissue": dict(tissue.most_common()), "target": dict(target.most_common(300))}


def find_models(query, family=None, assay=None, limit=25):
    """Free-text discovery over the index -- match `query` against tissue / target / assay / accession,
    ranked (exact > word-boundary > substring), QC-passed first. Returns ranked candidates. For a loose
    ORGAN word that is not literally in a cell-line name (e.g. 'liver' vs 'HepG2'), use
    find_models_by_organ, which asks ENCODE's biosample ontology."""
    idx = _load()
    q = (query or "").strip().lower()
    scored = []
    for acc, rec in idx.items():
        m, fam = rec["_meta"], rec["_family"]
        if family and fam != family.lower():
            continue
        if assay and assay.lower() not in (m.get("assay", "") or "").lower():
            continue
        best = 0
        for val in ((m.get("tissue", "") or "").lower(), (m.get("target", "") or "").lower(),
                    (m.get("assay", "") or "").lower(), acc.lower()):
            if not val or not q:
                continue
            if val == q:
                best = max(best, 100)
            elif val.startswith(q) or ((" " + q) in (" " + val)):
                best = max(best, 60)
            elif q in val:
                best = max(best, 30)
        if best:
            scored.append((best, 0 if (m.get("qc") == "passed") else 1,
                           {"accession": acc, "family": fam, "assay": m.get("assay"),
                            "target": m.get("target"), "tissue": m.get("tissue"), "qc": m.get("qc")}))
    scored.sort(key=lambda x: (-x[0], x[1]))
    return [r for _, _, r in scored[:limit]]


def find_models_by_organ(organ, family=None, assay_title=None, limit=25, timeout=30):
    """Ontology-aware discovery: ask ENCODE for experiments in an organ (biosample `organ_slims`),
    then keep the ones we have MODELS for. Catches 'liver' -> HepG2 / hepatocyte / liver tissue that a
    tissue-string match on our table would miss. Returns our matching models (accession, family, ...)."""
    import json as _json
    import urllib.parse
    import urllib.request
    q = [("type", "Experiment"), ("format", "json"), ("limit", "300"),
         ("biosample_ontology.organ_slims", organ)]
    if assay_title:
        q.append(("assay_title", assay_title))
    url = "https://www.encodeproject.org/search/?" + urllib.parse.urlencode(q)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        d = _json.load(r)
    idx = _load()
    out = []
    for h in d.get("@graph", []) or []:
        acc = h.get("accession")
        rec = idx.get(acc)
        if not rec:
            continue
        if family and rec["_family"] != family.lower():
            continue
        m = rec["_meta"]
        out.append({"accession": acc, "family": rec["_family"], "assay": m.get("assay"),
                    "tissue": m.get("tissue"), "target": m.get("target"), "qc": m.get("qc"),
                    "encode_biosample": h.get("biosample_summary")})
        if len(out) >= limit:
            break
    return out


# The slim facets ENCODE exposes on a biosample -- widened past organ so "immune B cell" (cell/system),
# "fetal brain" (developmental), etc. resolve, not just organs.
_SLIM_TYPES = ["organ_slims", "cell_slims", "system_slims", "developmental_slims"]

# Deterministic, auditable free-text -> slim-vocabulary map (the L1 query resolver; no ML). Maps loose
# words that are NOT themselves ENCODE slim values onto ones that are. Extend as real queries miss.
_SLIM_SYNONYMS = {
    "hepatic": "liver", "hcc": "liver", "hepatocellular": "liver", "hepatocyte": "liver",
    "bile organ": "liver", "b cell": "B cell", "b-cell": "B cell", "lymphoblast": "B cell",
    "lymphoblastoid": "B cell", "immune": "immune system", "immune cell": "immune system",
    "cardiac": "heart", "cardiomyocyte": "heart", "myocardial": "heart",
    "renal": "kidney", "neural": "brain", "neuronal": "brain", "neuron": "brain",
    "cortical": "brain", "cerebral": "brain", "pulmonary": "lung", "mammary": "breast",
    "gut": "intestine", "bowel": "intestine", "colonic": "large intestine",
}
# loose life-stage words -> ENCODE life_stage vocabulary
_LIFESTAGE_SYNONYMS = {"fetal": "embryonic", "foetal": "embryonic", "prenatal": "embryonic"}


def _norm_terms(term):
    """Candidate slim values for a free-text term: the term itself + synonym-map expansions (exact,
    then substring). De-duplicated, order-preserving."""
    if not term:
        return []
    t = term.strip().lower()
    cands = [term.strip()]
    if t in _SLIM_SYNONYMS:
        cands.append(_SLIM_SYNONYMS[t])
    for k, v in _SLIM_SYNONYMS.items():
        if k in t:
            cands.append(v)
    seen, out = set(), []
    for c in cands:
        if c.lower() not in seen:
            seen.add(c.lower()); out.append(c)
    return out


def find_models_by_biosample(term, family=None, assay_title=None, organism=None,
                             life_stage=None, slim_types=None, limit=25, timeout=30):
    """Discover models from a free-text tissue/cell/system/developmental description, WITHOUT the
    accession -- the L1 discovery path. Widens find_models_by_organ past organ_slims to cell/system/
    developmental slims, adds optional organism + life_stage facets, and resolves loose words via a
    deterministic synonym map. The ONTOLOGY decides membership (an ENCODE slim-facet set-intersection);
    similarity is never involved. Every returned model carries `authorized_by` = the exact facet+value
    that admitted it, so the choice is auditable. E.g. 'immune B cell' -> GM12878 via cell_slims;
    'mouse liver' -> mouse liver models via organism; 'hepatic' -> HepG2 via synonym(liver)+organ_slims.

    Returns [{accession, family, assay, tissue, target, qc, authorized_by:{facet,value},
    encode_biosample}]. Empty list = honest miss (no slim matched); the caller should abstain, never
    fall back to a fuzzy guess."""
    import json as _json
    import urllib.parse
    import urllib.request
    slim_types = slim_types or _SLIM_TYPES
    ls = life_stage
    if ls and ls.strip().lower() in _LIFESTAGE_SYNONYMS:
        ls = _LIFESTAGE_SYNONYMS[ls.strip().lower()]
    # The ENCODE BPNet-family Atlas index is human-only, and the tables carry no per-model organism
    # column -- so a non-human request is an honest EMPTY, never a silent human fallback. (When
    # non-human models are added, replace this with a real per-model organism filter.)
    if organism and organism.strip().lower() not in ("homo sapiens", "human", "h. sapiens", "hsapiens"):
        return []
    idx = _load()
    out, seen = [], set()
    for cand in _norm_terms(term):
        before = len(out)
        for st in slim_types:
            # Only the slim facet (+ assay_title) go in the ENCODE query -- both are known-good facets.
            # organism / life_stage are filtered CLIENT-SIDE from biosample_summary (their deep-path
            # facets 404 the whole query, which would silently zero the result).
            q = [("type", "Experiment"), ("format", "json"), ("limit", "300"),
                 ("biosample_ontology." + st, cand)]
            if assay_title:
                q.append(("assay_title", assay_title))
            url = "https://www.encodeproject.org/search/?" + urllib.parse.urlencode(q)
            try:
                req = urllib.request.Request(url, headers={"Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    d = _json.load(r)
            except Exception:
                continue   # invalid slim value for this facet -> ENCODE 404s; just try the next
            for h in d.get("@graph", []) or []:
                acc = h.get("accession")
                rec = idx.get(acc)
                if not rec or acc in seen:
                    continue
                if family and rec["_family"] != family.lower():
                    continue
                bs = h.get("biosample_summary") or ""
                if ls and bs and ls.lower() not in bs.lower():
                    continue   # best-effort life-stage filter (only when the summary carries it)
                m = rec["_meta"]
                seen.add(acc)
                out.append({"accession": acc, "family": rec["_family"], "assay": m.get("assay"),
                            "tissue": m.get("tissue"), "target": m.get("target"), "qc": m.get("qc"),
                            "authorized_by": {"facet": "biosample_ontology." + st, "value": cand},
                            "encode_biosample": h.get("biosample_summary")})
                if len(out) >= limit:
                    return out
            if len(out) > before:
                break   # this candidate matched a slim -> don't spend calls on the other slim types
    return out
