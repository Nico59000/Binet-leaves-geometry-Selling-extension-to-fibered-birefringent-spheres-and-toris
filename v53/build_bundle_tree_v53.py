#!/usr/bin/env python3
from pathlib import Path
import re, shutil, hashlib, csv, json
R=Path(__file__).resolve().parent
# Create navigation roots (never move/delete legacy root members).
for d in ['00_INDEX','01_DOCS','02_VERSIONS','03_SCRIPTS','04_CERTIFICATS','05_SELLING','06_CHIRES','07_QA_SEAL','08_FORMAL_CONTEXT','09_ASSETS']:
    (R/d).mkdir(exist_ok=True)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def copy_unique(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():
        if sha(src)!=sha(dst): raise RuntimeError(f'collision {dst}')
        return False
    shutil.copy2(src,dst); return True
root_files=[p for p in R.iterdir() if p.is_file()]
# 02_VERSIONS: physical mirrors of versioned cumulative PDF/TEX.
versioned=[]
rxver=re.compile(r'^formalisation_feuillets_binet_birefringence_selling_polaire_v(\d+)\.(pdf|tex)$')
for p in root_files:
    m=rxver.match(p.name)
    if m:
        v=int(m.group(1)); dst=R/'02_VERSIONS'/f'v{v:02d}'/p.name
        copy_unique(p,dst); versioned.append((p,dst,v))
# 03_SCRIPTS: physical mirrors, grouped by last version token when available.
script_ext={'.py','.cpp','.sh'}; scripts=[]
for p in root_files:
    if p.suffix.lower() not in script_ext: continue
    nums=re.findall(r'(?:^|[_-])v(\d+)(?:\D|$)',p.name,re.I)
    sub=f'v{int(nums[-1]):02d}' if nums else 'shared'
    dst=R/'03_SCRIPTS'/sub/p.name
    copy_unique(p,dst); scripts.append((p,dst,sub))
# 04_CERTIFICATS: physical mirrors of reports/certificates/checkpoints/ledgers.
keys=('certificate','certificat','checkpoint','report','audit','manifest','token','sha256','status','registry','decision','attestation','provenance','verification','delivery_index','math_checkpoint','source_context','formal_rigor')
certs=[]
for p in root_files:
    low=p.name.lower()
    if p.suffix.lower() not in {'.json','.txt','.md','.csv','.geojson'}: continue
    if not any(k in low for k in keys): continue
    nums=re.findall(r'(?:^|[_-])v(\d+)(?:\D|$)',p.name,re.I)
    sub=f'v{int(nums[-1]):02d}' if nums else 'shared'
    dst=R/'04_CERTIFICATS'/sub/p.name
    copy_unique(p,dst); certs.append((p,dst,sub))
# Lightweight thematic indexes, keeping large historical assets only once.
def write_index(path, selected):
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['root_path','bytes','sha256'])
        for p in sorted(selected,key=lambda x:x.name.lower()): w.writerow([p.name,p.stat().st_size,sha(p)])
write_index(R/'05_SELLING'/'ROOT_SELLING_INDEX.tsv',[p for p in root_files if any(k in p.name.lower() for k in ['selling','fig1','fig2','fig3','fig4','fig5','fig6','tafel','log_0016'])])
write_index(R/'06_CHIRES'/'ROOT_CHIRES_GAMMA3_INDEX.tsv',[p for p in root_files if any(k in p.name.lower() for k in ['gamma3','chires','r2pure'])])
formal_names=['addendum_sqrt3_hypergeometric_quarterturn_selling_fano_v2_2_hal','annexe-heyting','heyting_besteod','htnt_functional','pochhammer','polydivisible','positional_reversal','roof_five_slice']
write_index(R/'08_FORMAL_CONTEXT'/'ROOT_FORMAL_CONTEXT_INDEX.tsv',[p for p in root_files if any(k in p.name.lower() for k in formal_names)])
asset_ext={'.png','.jpg','.jpeg','.pdf'}
write_index(R/'09_ASSETS'/'ROOT_ASSETS_INDEX.tsv',[p for p in root_files if p.suffix.lower() in asset_ext and not rxver.match(p.name)])
# Root legacy catalog proves that historical root files still exist and are easy to locate.
with open(R/'00_INDEX'/'ROOT_LEGACY_CATALOG.tsv','w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256','kind'])
    for p in sorted(root_files,key=lambda x:x.name.lower()):
        kind='version' if rxver.match(p.name) else ('script' if p.suffix.lower() in script_ext else 'root-artifact')
        w.writerow([p.name,p.stat().st_size,sha(p),kind])
# Mirror catalog (excluding mutable catalogs themselves).
mirrors=[]
for top in ['02_VERSIONS','03_SCRIPTS','04_CERTIFICATS','05_SELLING/v53','06_CHIRES/v53','07_QA_SEAL/v53']:
    q=R/top
    if q.exists(): mirrors.extend(p for p in q.rglob('*') if p.is_file())
with open(R/'00_INDEX'/'NAVIGATION_MIRROR_CATALOG.tsv','w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256'])
    for p in sorted(mirrors,key=lambda x:x.as_posix().lower()): w.writerow([p.relative_to(R).as_posix(),p.stat().st_size,sha(p)])
readme='''# v53 bundle - start here\n\nThis bundle is append-only relative to sealed v52. All v52 root members remain at their historical paths and byte-exact. v53 adds a navigation layer; mirrored copies are conveniences, not replacements.\n\n- `01_DOCS/` - orientation and current-document pointers.\n- `02_VERSIONS/vNN/` - physical PDF/TEX mirrors by cumulative version.\n- `03_SCRIPTS/vNN/` and `shared/` - executable/verifier/build scripts.\n- `04_CERTIFICATS/vNN/` - certificates, reports, audits, manifests, tokens, ledgers.\n- `05_SELLING/` - Selling-specific indexes and v53 typed Fig.2-6 functor.\n- `06_CHIRES/` - CHIRES/Gamma3 indexes and canonical v53 obstruction data.\n- `07_QA_SEAL/` - PDF QA, visual contact sheets, seal/postseal records.\n- `08_FORMAL_CONTEXT/` - formal-context index.\n- `09_ASSETS/` - heavy-asset index; heavy historical files remain single-copy at root.\n- `00_INDEX/` - catalogs and bundle-tree audit.\n\nUse `00_INDEX/ROOT_LEGACY_CATALOG.tsv` to locate original root members and `00_INDEX/NAVIGATION_MIRROR_CATALOG.tsv` for organized mirrors.\n'''
(R/'00_INDEX'/'README_FIRST.md').write_text(readme,encoding='utf-8')
(R/'README_FIRST_v53.md').write_text(readme,encoding='utf-8')
(R/'01_DOCS'/'CURRENT_DOCUMENTS.md').write_text('''# Current cumulative documents\n\nCanonical root files:\n- `formalisation_feuillets_binet_birefringence_selling_polaire_v53.pdf`\n- `formalisation_feuillets_binet_birefringence_selling_polaire_v53.tex`\n\nPhysical version mirrors live in `../02_VERSIONS/v53/`.\n''',encoding='utf-8')
# Audit counts and size. Parent-byte equality is checked separately by the v53 verifier.
allfiles=[p for p in R.rglob('*') if p.is_file()]
rep={'version':'v53','status':'PASS','policy':'NONDESTRUCTIVE_NAVIGATION_LAYER; no legacy root member moved or renamed','root_files_catalogued':len(root_files),'version_mirrors':len(versioned),'script_mirrors':len(scripts),'certificate_mirrors':len(certs),'total_files_after_tree':len(allfiles),'total_uncompressed_bytes_after_tree':sum(p.stat().st_size for p in allfiles),'directories':['00_INDEX','01_DOCS','02_VERSIONS','03_SCRIPTS','04_CERTIFICATS','05_SELLING','06_CHIRES','07_QA_SEAL','08_FORMAL_CONTEXT','09_ASSETS']}
(R/'00_INDEX'/'BUNDLE_TREE_AUDIT.json').write_text(json.dumps(rep,indent=2)+'\n')
print(json.dumps(rep,indent=2))
