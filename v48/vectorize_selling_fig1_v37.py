#!/usr/bin/env python3
import csv, json, math, hashlib, subprocess, sys
from pathlib import Path
from collections import defaultdict
import cv2, numpy as np
from skimage.morphology import skeletonize
from scipy.ndimage import convolve

OUT=Path(__file__).resolve().parent
SRC=OUT/'selling_fig1_full_600.png'
PDF=OUT/'LOG_0016.pdf'
EXPECTED_PDF='423da98bd897c5dcb2b3556f4c9dc1367a4f12263d247608cfb2ba7a00f2e47b'
EXPECTED_EMBED='9116d054ae3267ae5f8a781135e434db659c4c8a6ff162cc5b6b2a12eee4c437'
EXPECTED_CROP='ff8583eb4372d9c4c1a285010e925173c78281c1cccf62888bdae47bca2261cf'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(PDF)==EXPECTED_PDF, 'LOG_0016.pdf hash mismatch'
assert sha(SRC)==EXPECTED_CROP, 'canonical crop hash mismatch'
# Verify embedded image bytes independently via pdfimages.
tmp=OUT/'_v37_pdfimages'; tmp.mkdir(exist_ok=True)
for q in tmp.glob('*'): q.unlink()
p=subprocess.run(['pdfimages','-all',str(PDF),str(tmp/'img')],capture_output=True,text=True)
assert p.returncode==0
imgs=list(tmp.glob('img-*'))
assert len(imgs)==1 and sha(imgs[0])==EXPECTED_EMBED

im600=cv2.imread(str(SRC),cv2.IMREAD_GRAYSCALE)
assert im600.shape==(5400,4450)
im=cv2.resize(im600,(2225,2700),interpolation=cv2.INTER_AREA)
cv2.imwrite(str(OUT/'selling_fig1_gray_300_v37.png'),im)
H,W=im.shape

# Long-source-stroke mask. Compact typography is excluded unless it belongs to a
# long/large connected component; this deliberately favours precision over recall.
b=(im<205).astype(np.uint8)
num,lab,stats,cent=cv2.connectedComponentsWithStats(b,8)
keep=np.zeros_like(b)
for i in range(1,num):
    x,y,w,h,a=map(int,stats[i]); mx=max(w,h); mn=min(w,h)
    elongated=mx>=75 and a>=35 and (mx/(mn+1)>=1.6 or a>=250)
    if elongated: keep[lab==i]=1
keep[:,:110]=0; keep[:180,:900]=0; keep[2550:,:]=0
keep=cv2.morphologyEx(keep,cv2.MORPH_CLOSE,np.ones((3,3),np.uint8),iterations=1)
cv2.imwrite(str(OUT/'selling_fig1_longstroke_mask_v37.png'),np.where(keep,0,255).astype(np.uint8))

# Skeleton graph.
sk=skeletonize(keep.astype(bool))
ker=np.ones((3,3),np.uint8); ker[1,1]=0
degpix=convolve(sk.astype(np.uint8),ker,mode='constant',cval=0)
nodepix=sk & (degpix!=2)
npix=cv2.dilate(nodepix.astype(np.uint8),np.ones((3,3),np.uint8),iterations=1)
nn,nlab,nstats,ncent=cv2.connectedComponentsWithStats(npix,8)
nodes_raw=[]
for i in range(1,nn):
    x,y,w,h,a=map(int,nstats[i]); cx,cy=ncent[i]
    nodes_raw.append({'raw_id':i-1,'x300':float(cx),'y300':float(cy),'bbox':[x,y,w,h],'area':a})
node_id_r=np.zeros_like(nlab,np.int32)
for compact,i in enumerate(range(1,nn),start=1): node_id_r[nlab==i]=compact
edgepix=sk & (node_id_r==0)
ne,elab,estats,ecent=cv2.connectedComponentsWithStats(edgepix.astype(np.uint8),8)
nbr=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
def order_component(coords):
    S=set(map(tuple,coords.tolist()))
    adj={p:[(p[0]+dy,p[1]+dx) for dy,dx in nbr if (p[0]+dy,p[1]+dx) in S] for p in S}
    ends=[p for p,v in adj.items() if len(v)<=1]
    start=ends[0] if ends else next(iter(S)); path=[start]; prev=None; cur=start; seen={cur}
    while True:
        cand=[q for q in adj[cur] if q!=prev and q not in seen]
        if not cand: break
        q=cand[0]; path.append(q); seen.add(q); prev,cur=cur,q
    return path

# stroke width field for classification
DT=cv2.distanceTransform(b,cv2.DIST_L2,5)
raw_edges=[]
for i in range(1,ne):
    x,y,w,h,a=map(int,estats[i]);
    if a<8: continue
    ys,xs=np.where(elab==i); coords=np.column_stack([ys,xs])
    comp=(elab==i).astype(np.uint8)
    ring=cv2.dilate(comp,np.ones((3,3),np.uint8),1).astype(bool) & ~comp.astype(bool)
    ids=sorted(set(int(v) for v in node_id_r[ring] if v>0))
    if len(ids)!=2: continue
    path=order_component(coords)
    xy=np.array([[p[1],p[0]] for p in path],dtype=np.float32)
    if len(xy)<2: continue
    plen=float(np.sum(np.hypot(np.diff(xy[:,0]),np.diff(xy[:,1]))))
    if plen<20: continue
    eps=max(1.2,min(4.0,plen/100))
    approx=cv2.approxPolyDP(xy.reshape(-1,1,2),eps,False).reshape(-1,2)
    chord=float(np.linalg.norm(approx[-1]-approx[0])) if len(approx)>1 else 0.0
    tort=plen/max(chord,1e-9)
    sam=[]
    for pxy in xy[::max(1,len(xy)//150)]:
        xx=int(round(float(pxy[0]))); yy=int(round(float(pxy[1])))
        if 2<=xx<W-2 and 2<=yy<H-2: sam.append(float(np.max(DT[yy-1:yy+2,xx-1:xx+2])))
    wmed=float(np.median(sam)) if sam else 0.0
    raw_edges.append({'raw_id':i-1,'raw_nodes':[j-1 for j in ids],'length_px300':plen,'chord_px300':chord,'tortuosity':tort,'stroke_radius_med_px300':wmed,'polyline300':[[float(px),float(py)] for px,py in approx]})

# Select source-anchored vector segments. This is deliberately not an exhaustive
# claim: weak/dotted layers are handled separately as candidates.
solid=[]; serp=[]
for e in raw_edges:
    if e['length_px300']<30: continue
    if e['stroke_radius_med_px300']>=2.7 and e['length_px300']>=50:
        layer='reinforced'; status='AUTO_HIGH_CONFIDENCE_SOURCE_STROKE'
    elif e['tortuosity']>=1.18 and e['stroke_radius_med_px300']<2.7:
        layer='serpentine'; status='AUTO_CANDIDATE_CURVED/WAVY_SOURCE_STROKE'
    else:
        layer='ordinary'; status='AUTO_HIGH_CONFIDENCE_SOURCE_STROKE'
    ee=dict(e); ee['layer']=layer; ee['status']=status
    (serp if layer=='serpentine' else solid).append(ee)

# Compact node set induced by ordinary+reinforced source strokes.
used=defaultdict(list)
for e in solid:
    for rn in e['raw_nodes']: used[rn].append(e)
node_map={}
node_rows=[]
for new_id,rn in enumerate(sorted(used)):
    r=nodes_raw[rn]; node_map[rn]=f'N{new_id:03d}'
    layers=defaultdict(int)
    for e in used[rn]: layers[e['layer']]+=1
    node_rows.append({
        'node_id':node_map[rn], 'x_px300':round(r['x300'],3),'y_px300':round(r['y300'],3),
        'x_px600':round(2*r['x300'],3),'y_px600':round(2*r['y300'],3),
        'degree_solid':len(used[rn]),'ordinary_degree':layers['ordinary'],'reinforced_degree':layers['reinforced'],
        'status':'SOURCE_STROKE_NODE_AUTO_CONFIRMED' if len(used[rn])>=2 else 'SOURCE_STROKE_ENDPOINT_AUTO_CONFIRMED'
    })

segment_rows=[]
def addseg(e, sid):
    p=e['polyline300']; n0=node_map.get(e['raw_nodes'][0],''); n1=node_map.get(e['raw_nodes'][1],'')
    return {
      'segment_id':sid,'layer':e['layer'],'node_u':n0,'node_v':n1,
      'length_px300':round(e['length_px300'],3),'tortuosity':round(e['tortuosity'],6),
      'stroke_radius_med_px300':round(e['stroke_radius_med_px300'],3),
      'polyline_px300':json.dumps([[round(x,2),round(y,2)] for x,y in p],separators=(',',':')),
      'status':e['status']
    }
for j,e in enumerate(solid): segment_rows.append(addseg(e,f'S{j:03d}'))
for j,e in enumerate(serp): segment_rows.append(addseg(e,f'W{j:03d}'))

# Symmetry/dotted candidate vectors from collinear groups of short LSD segments.
lsd=cv2.createLineSegmentDetector(cv2.LSD_REFINE_STD); L=lsd.detect(im)[0]
short=[]
for z in (L[:,0,:] if L is not None else []):
    x1,y1,x2,y2=map(float,z); le=math.hypot(x2-x1,y2-y1)
    if le<7: continue
    mx=(x1+x2)/2; my=(y1+y2)/2
    if mx<110 or my>2520 or (my<220 and mx<1300): continue
    ang=math.atan2(y2-y1,x2-x1); ang=ang+math.pi if ang<0 else ang
    dx,dy=math.cos(ang),math.sin(ang); nx,ny=-dy,dx
    rho=nx*mx+ny*my; t1=dx*x1+dy*y1; t2=dx*x2+dy*y2
    if t1>t2:t1,t2=t2,t1
    short.append((x1,y1,x2,y2,le,ang,rho,t1,t2))
# coarse grouping; only export strongest non-overlapping candidates, never call them confirmed axes
bins=defaultdict(list)
for s in short:
    bins[(round(math.degrees(s[5])/2.5),round(s[6]/7))].append(s)
cands=[]
for key,arr in bins.items():
    if len(arr)<4: continue
    # use bin centre geometry
    a=np.median([s[5] for s in arr]); dx,dy=math.cos(a),math.sin(a); nx,ny=-dy,dx
    rho=float(np.median([s[6] for s in arr])); ints=sorted((s[7],s[8]) for s in arr)
    merged=[]
    for lo,hi in ints:
        if not merged or lo>merged[-1][1]+3: merged.append([lo,hi])
        else: merged[-1][1]=max(merged[-1][1],hi)
    if len(merged)<4: continue
    tmin=min(x[0] for x in merged); tmax=max(x[1] for x in merged); span=tmax-tmin
    cov=sum(hi-lo for lo,hi in merged)/max(span,1)
    if span<350 or not(0.06<=cov<=0.30): continue
    x1=rho*nx+tmin*dx; y1=rho*ny+tmin*dy; x2=rho*nx+tmax*dx; y2=rho*ny+tmax*dy
    cands.append({'x1':x1,'y1':y1,'x2':x2,'y2':y2,'span':span,'coverage':cov,'pieces':len(merged),'angle_deg':math.degrees(a),'rho':rho})
# de-duplicate and take top 18 candidates; source asserts 6 axes but matching remains unpromoted
uniq=[]
for c in sorted(cands,key=lambda q:-(q['span']*q['pieces'])):
    ok=True
    for d in uniq:
        ad=abs(c['angle_deg']-d['angle_deg']); ad=min(ad,180-ad)
        if ad<2 and abs(c['rho']-d['rho'])<12: ok=False; break
    if ok: uniq.append(c)
    if len(uniq)>=18: break
for j,c in enumerate(uniq):
    segment_rows.append({
      'segment_id':f'Y{j:03d}','layer':'symmetry_candidate','node_u':'','node_v':'',
      'length_px300':round(c['span'],3),'tortuosity':1.0,'stroke_radius_med_px300':'',
      'polyline_px300':json.dumps([[round(c['x1'],2),round(c['y1'],2)],[round(c['x2'],2),round(c['y2'],2)]],separators=(',',':')),
      'status':'AUTO_DASHED/COLLINEAR_CANDIDATE_NOT_MATCHED_TO_SIX_SOURCE_AXES'
    })

# Stable closed-region candidates from the source ink. These are NOT promoted to
# historical fields; stability across morphology thresholds is recorded.
def regions(k):
    base=(im<200).astype(np.uint8)
    m=cv2.morphologyEx(base,cv2.MORPH_CLOSE,np.ones((k,k),np.uint8),iterations=1)
    m=cv2.dilate(m,np.ones((3,3),np.uint8),1); bg=1-m
    n,la,st,ce=cv2.connectedComponentsWithStats(bg,4); out=[]
    for i in range(1,n):
        x,y,w,h,a=map(int,st[i]); cx,cy=ce[i]
        touch=(x==0 or y==0 or x+w>=W or y+h>=H)
        if not touch and a>=1000 and w>=20 and h>=20 and 100<cx<2100 and 100<cy<2550:
            out.append({'label':i,'area':a,'cx':float(cx),'cy':float(cy),'w':w,'h':h,'lab':la})
    return out,m
r15,_=regions(15); r17,m17=regions(17); r19,_=regions(19)
# solid vector mask for boundary support score
solidmask=np.zeros_like(im,np.uint8)
for e in solid:
    pts=np.round(np.array(e['polyline300'])).astype(np.int32).reshape(-1,1,2)
    cv2.polylines(solidmask,[pts],False,1,5,cv2.LINE_AA)
solidD=cv2.dilate(solidmask,np.ones((9,9),np.uint8),1)
def match_stable(r, others):
    n=1
    for O in others:
        best=min(O,key=lambda q:math.hypot(q['cx']-r['cx'],q['cy']-r['cy']),default=None)
        if best and math.hypot(best['cx']-r['cx'],best['cy']-r['cy'])<=18 and 0.5<=best['area']/r['area']<=2.0: n+=1
    return n
field_rows=[]
for j,r in enumerate(r17):
    reg=(r['lab']==r['label']).astype(np.uint8); er=cv2.erode(reg,np.ones((5,5),np.uint8),1); bd=reg-er
    ov=float((solidD[bd.astype(bool)]>0).mean()) if bd.sum() else 0.0
    stable=match_stable(r,[r15,r19])
    if stable==3 and r['area']>=8000 and ov>=0.35: cls='STABLE_CLOSED_REGION_HIGH'
    elif stable==3 and r['area']>=2500 and ov>=0.30: cls='STABLE_CLOSED_REGION_MEDIUM'
    else: cls='MORPHOLOGICAL_CANDIDATE'
    field_rows.append({
      'region_id':f'R{j:03d}','centroid_x_px300':round(r['cx'],3),'centroid_y_px300':round(r['cy'],3),
      'area_px2_300':r['area'],'bbox_w_px300':r['w'],'bbox_h_px300':r['h'],
      'stability_k15_k17_k19':stable,'solid_boundary_support':round(ov,6),
      'classification':cls,'historical_field_status':'NOT_PROMOTED_WITHOUT_LAYER/INCIDENCE_MANUAL_CONFIRMATION'
    })

# CSV exports.
node_path=OUT/'selling_fig1_nodes_v37.csv'; seg_path=OUT/'selling_fig1_segments_v37.csv'; field_path=OUT/'selling_fig1_fields_v37.csv'
with node_path.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(node_rows[0].keys()) if node_rows else ['node_id']); w.writeheader(); w.writerows(node_rows)
with seg_path.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(segment_rows[0].keys()) if segment_rows else ['segment_id']); w.writeheader(); w.writerows(segment_rows)
with field_path.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(field_rows[0].keys()) if field_rows else ['region_id']); w.writeheader(); w.writerows(field_rows)

# Layer GeoJSON in image coordinates.
features=[]
for s in segment_rows:
    pts=json.loads(s['polyline_px300'])
    features.append({'type':'Feature','geometry':{'type':'LineString','coordinates':pts},'properties':{k:v for k,v in s.items() if k!='polyline_px300'}})
geo={'type':'FeatureCollection','crs_note':'pixel coordinates, canonical 300-dpi crop; origin top-left','features':features}
geo_path=OUT/'selling_fig1_layers_v37.geojson'; geo_path.write_text(json.dumps(geo,indent=2,ensure_ascii=False),encoding='utf-8')

# Diagnostic overlay.
rgb=cv2.cvtColor(im,cv2.COLOR_GRAY2BGR)
cols={'ordinary':(0,0,220),'reinforced':(0,120,255),'serpentine':(180,0,180),'symmetry_candidate':(220,100,0)}
for s in segment_rows:
    pts=np.round(np.array(json.loads(s['polyline_px300']))).astype(np.int32).reshape(-1,1,2)
    cv2.polylines(rgb,[pts],False,cols[s['layer']],2 if s['layer']=='reinforced' else 1,cv2.LINE_AA)
for n in node_rows:
    cv2.circle(rgb,(round(n['x_px300']),round(n['y_px300'])),2,(0,180,0),-1)
cv2.imwrite(str(OUT/'selling_fig1_layers_overlay_v37.png'),rgb)

manifest={
 'source_pdf':PDF.name,'source_pdf_sha256':sha(PDF),'embedded_image_sha256':EXPECTED_EMBED,
 'canonical_crop':SRC.name,'canonical_crop_sha256':sha(SRC),'crop_dimensions_600':[4450,5400],
 'tables':{
  node_path.name:{'rows':len(node_rows),'sha256':sha(node_path)},
  seg_path.name:{'rows':len(segment_rows),'sha256':sha(seg_path)},
  field_path.name:{'rows':len(field_rows),'sha256':sha(field_path)},
  geo_path.name:{'features':len(features),'sha256':sha(geo_path)}},
 'layer_counts':dict((k,sum(1 for s in segment_rows if s['layer']==k)) for k in ['ordinary','reinforced','serpentine','symmetry_candidate']),
 'field_region_class_counts':dict((k,sum(1 for r in field_rows if r['classification']==k)) for k in ['STABLE_CLOSED_REGION_HIGH','STABLE_CLOSED_REGION_MEDIUM','MORPHOLOGICAL_CANDIDATE']),
 'decision':'SOURCE_BYTES_PROVEN; VECTOR_LAYERS_SOURCE_ANCHORED; COMPLETE_HISTORICAL_FIELD_INCIDENCE_NOT_YET_PROMOTED',
 'guard':'automatic region closure and dashed-line candidates are diagnostics, not historical incidence facts until manually/source-geometrically matched.'
}
(OUT/'selling_fig1_vectorization_manifest_v37.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(manifest,indent=2,ensure_ascii=False))
