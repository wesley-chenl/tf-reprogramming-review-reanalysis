import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def tsv(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def write(p, rows, fields=None):
    p.parent.mkdir(parents=True,exist_ok=True)
    fields=fields if fields is not None else list(dict.fromkeys(k for r in rows for k in r))
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def json_write(p,x):
    p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def js(p):
    return json.loads(p.read_text(encoding="utf-8"))


def verify_lock(root):
    rows=[]
    for line in (root/"SCORES_LOCK.sha256").read_text(encoding="utf-8-sig").splitlines():
        digest,name=line.split(maxsplit=1)
        name=name.lstrip("*")
        p=root/name
        if not p.resolve().is_relative_to(root.resolve()) or sha(p)!=digest:
            raise ValueError(f"Score lock violation: {name}")
        rows.append(name)
    return rows
