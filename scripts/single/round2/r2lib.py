"""Shared paths and durable structured records for the three-target benchmark."""
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(os.environ.get("HAMMELMAN_ROUND2", Path(__file__).resolve().parents[1])).resolve()
BASE = Path(os.environ.get("HAMMELMAN_BASE", ROOT.parent / "base")).resolve()
FOLLOW = ROOT.parent / "qc"
AUDIT = ROOT.parent / "annotation-audit"
RUNTIME = ROOT / "toolchain"
BASE_ENV = Path(os.environ.get("HAMMELMAN_BASE_ENV", str(BASE / "environments/posix_runtime/bench")))
R_ENV = Path(os.environ.get("HAMMELMAN_R_LIB", str(RUNTIME / "Rlib")))
ANA = Path(os.environ.get("HAMMELMAN_ANANSE_ENV", str(RUNTIME / "ananse")))
SCRATCH = ROOT / "scratch"
TASKS = ("liver", "cardiomyocyte", "betapanc")
METHODS = ("DE", "TRIAGE", "VIPER", "Taiji", "ANANSE", "CompForce RRA", "IMCell")
CAP = {"liver": "Liver", "cardiomyocyte": "Cardiomyocytes", "betapanc": "Pancreatic β cells"}
RNA_SRC = ("GSM1866625", "GSM1866626", "GSM1866627")
RNA_TGT = {
    "liver": ("ENCSR096STK_bio1", "ENCSR096STK_bio2"),
    "cardiomyocyte": ("GSM2525582", "GSM2525583"),
    "betapanc": ("GSM2666438", "GSM2666440"),
}
ATAC_SRC = ("GSM1867737", "GSM1867738")
ATAC_TGT = {"liver": "ENCSR609OHJ_bio1", "cardiomyocyte": "GSM2525591", "betapanc": "GSM2666437"}


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sg():
    return dt.datetime.now(ZoneInfo("Asia/Singapore")).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while True:
            block = stream.read(8 * 1024 * 1024)
            if not block:
                return digest.hexdigest()
            digest.update(block)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.partial")
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, columns=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if columns is None:
        columns = list(dict.fromkeys(k for row in rows for k in row))
    temporary = path.with_name(path.name + f".{os.getpid()}.partial")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=columns, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def input_record(path, role, digest=None, **extra):
    path = Path(path)
    stat = path.stat()
    return {"path": str(path), "role": role, "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "sha256": digest or sha(path), **extra}


def linux_path(value):
    return Path(f"/mnt/{value[0].lower()}/{value[3:].replace(chr(92), '/')}" if len(value) > 2 and value[1] == ":" else value)


def run_env():
    env = os.environ.copy()
    env.update({key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                   "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")})
    env["PATH"] = str(BASE_ENV / "bin") + ":" + env.get("PATH", "")
    env["R_LIBS_USER"] = str(R_ENV)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["MPLCONFIGDIR"] = str(RUNTIME / "matplotlib")
    env["XDG_CACHE_HOME"] = str(RUNTIME / "cache")
    env["XDG_CONFIG_HOME"] = str(RUNTIME / "config_home")
    return env
