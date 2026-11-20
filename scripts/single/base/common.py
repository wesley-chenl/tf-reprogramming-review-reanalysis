from pathlib import Path
import os,sys,json,csv,hashlib,datetime,urllib.request,time
OUT=Path(os.environ.get('HAMMELMAN_BASE', Path(__file__).resolve().parents[1])).resolve()
ROOT=OUT.parent
AUDIT=ROOT/'gemini_independent_audit_20261117_002158'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='4'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
_FROZEN_HASHES=None
def digest(p):
 global _FROZEN_HASHES
 p=Path(p)
 if _FROZEN_HASHES is None:
  cp=OUT/'config/FROZEN_INPUT_HASHES.json'
  _FROZEN_HASHES=json.loads(cp.read_text()) if cp.exists() else {}
 rec=_FROZEN_HASHES.get(str(p.resolve()))
 if rec:
  st=p.stat()
  if st.st_size!=rec['bytes'] or st.st_mtime_ns!=rec['mtime_ns']:raise RuntimeError('Frozen input changed: '+str(p))
  return rec['sha256']
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def dump(name,data):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def tsv(name,rows,fields=None):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True)
 if fields is None:fields=list(dict.fromkeys(k for r in rows for k in r))
 with open(p,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def fetch(url,relative,timeout=35,max_bytes=None):
 p=OUT/relative;p.parent.mkdir(parents=True,exist_ok=True);r={'url':url,'path':str(p),'started_UTC':now()}
 try:
  if not p.exists():
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Hammelman focused computational reanalysis'}),timeout=timeout) as x,open(str(p)+'.partial','wb') as f:
    r.update(http_status=x.status,headers=dict(x.headers));n=0
    while b:=x.read(1024*1024):
     n+=len(b)
     if max_bytes and n>max_bytes:raise ValueError('declared download-size limit')
     f.write(b)
   Path(str(p)+'.partial').replace(p)
  r.update(status='downloaded',bytes=p.stat().st_size,sha256=digest(p))
 except Exception as e:r.update(status='failed',error=repr(e))
 r['ended_UTC']=now();print(json.dumps({k:v for k,v in r.items() if k!='headers'}),flush=True);return r
def progress(message,phase='running'):
 dump('RUN_STATUS.json',{'phase':phase,'updated_UTC':now(),'message':message})
 p=OUT/'PROGRESS.md'
 with open(p,'a',encoding='utf-8') as f:f.write(f'\n- {now()}: {message}\n')
