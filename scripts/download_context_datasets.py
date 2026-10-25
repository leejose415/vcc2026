#!/usr/bin/env python3
"""Download a recorded source manifest with resumable transfers and checksum validation."""
import concurrent.futures, hashlib, json, pathlib, shutil, subprocess, sys, time
ROOT=pathlib.Path(__file__).resolve().parents[1]
manifest=pathlib.Path(sys.argv[1])
rows=json.loads(manifest.read_text())
def download(r):
 p=ROOT/r['path'];p.parent.mkdir(parents=True,exist_ok=True)
 part=p.with_name(p.name+'.part')
 if part.with_name(part.name+'.aria2').exists():raise RuntimeError('Parallel transfer state exists; resume with aria2, not curl: '+str(part))
 def valid(f):
  if not f.exists() or f.stat().st_size!=r['size']:return False
  h=hashlib.new(r.get('algorithm','md5'))
  with f.open('rb') as stream:
   for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
  return h.hexdigest()==r['checksum']
 if valid(p):print('VERIFIED existing',r['path'],flush=True);return
 if p.exists():raise RuntimeError('Existing file does not match manifest: '+str(p))
 if shutil.disk_usage(ROOT).free < r['size']-(part.stat().st_size if part.exists() else 0)+10*2**30:raise RuntimeError('Insufficient space with 10 GiB reserve')
 print('START',r['path'],r['size'],flush=True)
 subprocess.run(['curl','--fail','--location','--silent','--show-error','--retry','5','--retry-delay','5','--connect-timeout','30','--speed-limit','1024','--speed-time','120','--continue-at','-','--output',str(part),r['url']],check=True)
 if not valid(part):raise RuntimeError('Checksum or size mismatch: '+str(part))
 part.rename(p);print('VERIFIED',r['path'],flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 futures=[pool.submit(download,r) for r in rows]
 errors=[]
 for f in concurrent.futures.as_completed(futures):
  try:f.result()
  except Exception as e:errors.append(str(e));print('ERROR',e,flush=True)
print('FINISHED',len(rows)-len(errors),'verified;',len(errors),'errors',flush=True)
if errors:sys.exit(1)
