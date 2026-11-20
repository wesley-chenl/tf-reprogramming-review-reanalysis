import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def record(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def table(path, rows, fields=None):
    fields = fields if fields is not None else list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def check_lock(root, name):
    lock = load(root / "locks" / (name + "_LOCK.json"))
    for relative, digest in lock["files"].items():
        if sha(root / relative) != digest:
            raise ValueError(f"Frozen {name} artifact changed: {relative}")
    return lock
